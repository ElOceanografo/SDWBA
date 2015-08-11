
using DataFrames, Distributions, PyPlot
import Base: copy, length

type Scatterer{T}
	r::Array{T, 2}
	a::Array{T, 1}
	h::Array{T, 1}
	g::Array{T, 1}
	f0::T
end
copy(s::Scatterer) = Scatterer(s.r, s.a, s.h, s.g, s.f0)

function rescale(s::Scatterer; scale=1.0, radius=1.0, x=1.0, y=1.0, z=1.0)
	s = copy(s)
	M = diagm([x, y, z]) * scale
	s.r = M * s.r
	s.a = s.a * scale * radius
	return s
end

length(s::Scatterer) = norm(s.r[:, 1] - s.r[:, end])

function rotate(s::Scatterer; roll=0.0, tilt=0.0, yaw=0.0)
	roll, tilt, yaw = deg2rad([roll, tilt, yaw])
	Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)]
	Ry = [cos(tilt) 0 sin(tilt); 0 1 0; -sin(tilt) 0 cos(tilt)]
	Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1]
	R = Rz * Ry * Rx
	s1 = copy(s)
	s1.r = R * s1.r
	return s1
end

function DWBAintegrand(s, rr, aa, gg, hh, k)
	r = vec(rr[:, 1] + s * diff(rr, 2))
	a = aa[1] + s * diff(aa)[1]
	g = gg[1] + s * diff(gg)[1]
	h = hh[1] + s * diff(hh)[1]
	alphatilt = acos(dot(k, r) / (norm(k) * norm(r)))
	betatilt = abs(alphatilt - pi/2)
	gammagamma = 1 / (g * h^2) + 1 / g - 2
	if abs(abs(betatilt) - pi/2) < 1e-10
		bessy = norm(k) * a / h
	else
		bessy = besselj1(2*norm(k) * a / h * cos(betatilt)) / cos(betatilt)
	end
	return (norm(k) / 4 * gammagamma * exp(2im * dot(k,r) / h) * a * 
		bessy * norm(diff(rr, 2)))
end

function form_function(s::Scatterer, k, phase_sd=0.0)
	fbs = 0 + 0im
	m, n = size(s.r)

	for i in 1:(n-1)
		f(x) = DWBAintegrand(x, s.r[:, i:i+1], s.a[i:i+1],
			s.g[i:i+1], s.h[i:i+1], k)
		fbs += quadgk(f, 0.0, 1.0)[1] * exp(im * randn() * phase_sd)
	end
	return fbs
end

function backscatter_xsection{T}(s::Scatterer{T}, k::Vector{T}, phase_sd=0.0)
	return abs2(form_function(s, k, phase_sd))
end

function target_strength{T}(s::Scatterer{T}, k::Vector{T}, phase_sd=0.0)
	return 10 * log10(backscatter_xsection(s, k, phase_sd))
end


###################################################################
krill_df = readtable("data/generic_krill_McGeehee1998.csv")
krill_df[:y] = 0.0

krill = Scatterer(
	convert(Array, krill_df[[:x, :y, :z]])',
	convert(Array, krill_df[:a]),
	convert(Array, krill_df[:h]),
	convert(Array, krill_df[:g]),
	120e3)

krill.r *= 1e-3
krill.a *= 1e-3
c = 1456.0
k = [0, 0, 2pi * -120e3 / c]


angles = [-90:270.0]
TS = zeros(angles)
TS_stochastic = zeros(angles)
for i in 1:length(angles)
	println(i)
	TS[i] = target_strength(rotate(krill, tilt=angles[i]), k)
	t = Float64[target_strength(rotate(krill, tilt=angles[i]), k, 0.7071)
		for j in 1:100]
	TS_stochastic[i] = mean(t)
end
plot(angles, TS)
plot(angles, TS_stochastic)


freqs = [10:1000] * 1e3
TSf = zeros(freqs)
for i in 1:length(freqs)
	k = [0, 0, -2pi * freqs[i] / c]
	TSf[i] = target_strength(krill, k)
end
semilogx(freqs/1e3, TSf)


nrep = 5000
TS120 = zeros(nrep)
TS153 = zeros(nrep)
yaw = Uniform(0, 360)
tilt = Normal(11, 6)
adcp_angle = deg2rad(30)
k120 = [0, 0, 2pi * 120e3 / c]
k153 = [sin(adcp_angle), 0, cos(adcp_angle)] * 2pi * 153e3 / c
krill_small = rescale(krill, scale=0.8)

for i in 1:nrep
	TS120[i] = target_strength(rotate(krill_small, yaw=rand(yaw), tilt=rand(tilt)),
		k120, 0.7071)
	TS153[i] = target_strength(rotate(krill_small, yaw=rand(yaw), tilt=rand(tilt)),
		k153, 0.7071)
end

bins = [-120:-65]
subplot(211)
plt.hist(TS120, bins, normed=true)
xlim(-120, -65)
subplot(212)
plt.hist(TS153, bins, normed=true)
xlim(-120, -65)

println(mean(TS120), " ", mean(TS153))
println(std(TS120), " ", std(TS153))
