library(numDeriv)

sound.speed.water <- 1460
density.water <- 1027
freq <- 120e3

krill.shape <- read.csv("generic_krill.csv")
krill.shape <- krill.shape[order(krill.shape$x), ]
krill.shape$sound.speed <- sound.speed.water * krill.shape$h
krill.shape$density <- density.water * krill.shape$g
# convert form meters to mm
krill.shape[c("x", "y", "z", "a")] <- krill.shape[c("x", "y", "z", "a")] * 1e-3

dx <- diff(krill.shape$x)
dy <- diff(krill.shape$y)
dz <- diff(krill.shape$z)
krill.shape$s <- c(0, cumsum(sqrt(dx^2 + dy^2 + dz^2)))

x_fun <- splinefun(krill.shape$s, krill.shape$x)
y_fun <- splinefun(krill.shape$s, krill.shape$y)
z_fun <- splinefun(krill.shape$s, krill.shape$z)
a_fun <- splinefun(krill.shape$s, krill.shape$a)
g_fun <- splinefun(krill.shape$s, krill.shape$g)
h_fun <- splinefun(krill.shape$s, krill.shape$h)

r_fun <- function(s) c(x_fun(s), y_fun(s), z_fun(s))
local_tangent <- function(s) jacobian(r_fun, s)
dot <- function(x, y) sum(x * y)
vecnorm <- function(x) sqrt(sum(x^2))

tilt <- pi/2
k <- 2 * pi * freq * c(cos(tilt), sin(tilt), 0)



integrand <- function(s, k) {
  loc_tan <- local_tangent(s)
  beta <- acos(dot(k, loc_tan) / (vecnorm(k) * vecnorm(loc_tan)))
  beta <- abs(beta - pi/2)
  (g_fun(s) - h_fun(s)) * exp(2i * dot(k, r_fun(s))) * 
    besselJ(2 * vecnorm(k) * a_fun(s) * cos(beta), 1) / cos(beta)
}

integrand.real <- Vectorize(function(s, k) Re(integrand(s, k)), vectorize.args=c("s"))
integrand.imag <- Vectorize(function(s, k) Im(integrand(s, k)), vectorize.args=c("s"))

f.bs.real <- integrate(integrand.real, 0, max(krill.shape$s), k=k)$value
f.bs.imag <- integrate(integrand.imag, 0, max(krill.shape$s), k=k)$value
10 * log10(abs(f.bs.real + f.bs.imag)^2)

par(mfrow=c(2, 1))
plot(y ~ x, krill.shape, ty='b', ylim=c(-0.001, 0.005))
lines(krill.shape$x, krill.shape$y + krill.shape$a, lty=2)
lines(krill.shape$x, krill.shape$y - krill.shape$a, lty=2)

ss <- seq(0, max(krill.shape$s), length.out=500)
lines(x_fun(ss), y_fun(ss), col='red')
lines(x_fun(ss), y_fun(ss) + a_fun(ss), col='blue')
lines(x_fun(ss), y_fun(ss) - a_fun(ss), col='blue')

plot(x_fun(ss), abs(integrand.real(ss, k) + integrand.imag(ss, k))^2, ty='l')
