
using DataFrames, PyPlot

krill_df = readtable("generic_krill_McGeehee1998.csv")
krill_df[:y] = 0.0

type Scatterer{T}
	r::Array{T, 2}
	a::Array{T, 1}
	h::Array{T, 1}
	g::Array{T, 1}
	f0::T
end

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


function DWBAintegrand(s, rr, aa, gg, hh, k)
	r = s * diff(rr, 2)
	a = s * diff(aa, 2)
	g = s * diff(gg, 2)
	h = s * diff(hh, 2)
	alphatilt = acos(adot(k, r) / (norm(k) * norm(r)))
	betatilt = abs(alphatilt - pi/2)
	
end

function form_function(s::Scatterer)
	fbs = 0 + 0im
	m, n = size(s.r)

	for i in 1:(n-1)
		a1 = s.a[i]
		a2 = s.a[i+1]
		r1 = s.r[]
end
