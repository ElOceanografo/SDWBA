
## functions based on McGeehee et al. 1998
norm <- function(x) sqrt(sum(x^2))
dot <- function(x, y) sum(x * y)
rad2deg <- function(x) x * 180 / pi
deg2rad <- function(x) x * pi / 180


DWBAintegrand <- function(s, args) {
  rx <- s * (args$r2[1] - args$r1[1]) + args$r1[1]
  ry <- s * (args$r2[2] - args$r1[2]) + args$r1[2]
  rz <- s * (args$r2[3] - args$r1[3]) + args$r1[3]
  r <- c(rx, ry, rz)
  a <- s * (args$a2 - args$a1) + args$a1
  h <- s * (args$h2 - args$h1) + args$h1
  g <- s * (args$g2 - args$g1) + args$g1
  gamgam <- 1 / (g * h^2) + 1 / g - 2
  if (abs(args$betatilt) < 1e-10) {
    bessy <- norm(args$k1) * a / h
  } else {
#     x <- 2 * norm(args$k1) * a / h * cos(args$betatilt)
    bessy <- besselJ(2 * norm(args$k1) * a / h * cos(args$betatilt), 1) / cos(args$betatilt)
  }
  return( norm(args$k1) / 4 * gamgam * exp(2i * dot(args$k1, r) / h) * a * 
            bessy * norm(args$r2 - args$r1))
}
DWBAintegrand.Re <- Vectorize(function(s, args) Re(DWBAintegrand(s, args)), 
                              vectorize.args=c("s"))
DWBAintegrand.Im <- Vectorize(function(s, args) Im(DWBAintegrand(s, args)),
                              vectorize.args=c("s"))


backscatter.form <- function(x, y, z, a, g, h, k) {
  fbs <- 0 + 0i
  for (j in 1:(length(x) - 1)) {
#     print(paste("Segment", j))
    r1 <- c(x[j], y[j], z[j])
    r2 <- c(x[j + 1], y[j + 1], z[j + 1])
    alphatilt <- acos(dot(k, (r2 - r1)) / (norm(k) * norm(r2 - r1)))
    betatilt <- abs(alphatilt - pi/2)
    
    args <- list(r1=r1, r2=r2, a1=a[j], a2=a[j+1], g1=g[j], g2=g[j+1], 
                 h1=h[j], h2=h[j+1], betatilt=betatilt, k1=k)
    fbs  <- fbs + integrate(DWBAintegrand.Re, 0, 1, args=args)$value +
      integrate(DWBAintegrand.Im, 0, 1, args=args)$value * 1i
  }
  return(fbs)
}

k.magnitude <- 2 * pi * 120e3 / 1480
tilt <- 45 * pi / 180
k <- k.magnitude * c(cos(tilt), sin(tilt), 0)
krill.shape <- read.csv("generic_krill.csv")
krill.shape[c('x', 'y', 'z', 'a')] <- krill.shape[c('x', 'y', 'z', 'a')] * 0.001

f.bs <- with(krill.shape, backscatter.form(x, y, z, a, g, h, k))
10 * log10(abs(f.bs)^2)

theta <- 0:360
TS <- rep(0, length(theta))
for (i in 1:length(theta)) {
  k <- k.magnitude * c(cos(deg2rad(theta[i])), 0, sin(deg2rad(theta[i])))
  fbs <- with(krill.shape, backscatter.form(x, y, z, a, g, h, k))
  TS[i] <- 10 * log10(abs(fbs)^2)
}
plot(theta, TS, ty='l')


freqs <- seq(50e3, 500e3, 1e3)
theta <- 5
theta <- theta * pi / 180
TS <- rep(0, length(freqs))
for (i in 1:length(freqs)) {
  k <- c(sin(theta), 0, cos(theta)) * freqs[i] * 2 * pi / 1480
  fbs <- with(krill.shape, backscatter.form(x, y, z, a, g, h, k))
  TS[i] <- 10 * log10(abs(fbs)^2)
}
plot(freqs, TS, ty='l')

