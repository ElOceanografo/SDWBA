

segment.form.function <- function(r, a, angle, beta.tilt, segment.length, sound.speed, 
                                  density, freq, sound.speed.water, density.water) {
  gamma.rho <- (density - density.water) / density
  gamma.kappa <- density.water * sound.speed.water^2 / (density * sound.speed^2)
  k1 <- 2 * pi * freq / sound.speed.water
  k1.vec <- k1 * c(cos(angle), 0, sin(angle))
  k2 <- 2 * pi * freq / sound.speed
  i <- complex()
  integrand <- a * (gamma.kappa - gamma.rho) * exp(-2i * sum(r * k1.vec)) *
    besselJ(2 * k2 * a * abs(cos(beta.tilt)), 1) / cos(beta.tilt)
  # abs(cos(beta.tilt)) above is to make sure argument is positive...need to make sure that's ok
#   print(integrand)
  return(k1 / 4 * integrand * segment.length)
}


backscatter.xsection <- function(animal.body, freq, angle, sound.speed.water,
                                 density.water, phase.sd=0) {
  # for now, angle starts on x-axis and goes counterclockwise
  f.total <- 0i
  angle <- angle * pi / 180
  dx <- diff(animal.body$x)
  dy <- diff(animal.body$y)
  dz <- diff(animal.body$z)
  segment.length <- sqrt(dx^2 + dy^2 + dz^2)
  beta.tilt <- angle - atan2(dy, dx)
  
  for (j in 1:(nrow(animal.body) - 1)) {
    r <- animal.body[j, c("x", "y", "z")]
    a <- animal.body[j, "a"]
    sound.speed <- animal.body[j, "sound.speed"]
    density <- animal.body[j, "density"]    
    f.total <- f.total + segment.form.function(r, a, angle, beta.tilt[j], segment.length[j],
      sound.speed, density, freq, sound.speed.water, density.water)
    f.total <- f.total * exp(1i * rnorm(1, 0, phase.sd))
  }
  return(abs(f.total)^2)
}


sound.speed.water <- 1460
density.water <- 1027

krill.shape <- read.csv("generic_krill.csv")
krill.shape <- krill.shape[order(krill.shape$x), ]
krill.shape$sound.speed <- sound.speed.water * krill.shape$h
krill.shape$density <- density.water * krill.shape$g
# convert form meters to mm
krill.shape[c("x", "y", "z", "a")] <- krill.shape[c("x", "y", "z", "a")] * 1e-3

plot(y ~ x, krill.shape, ty='b', ylim=c(-0.001, 0.005))
lines(krill.shape$x, krill.shape$y + krill.shape$a, lty=2)
lines(krill.shape$x, krill.shape$y - krill.shape$a, lty=2)

sigma.bs <- backscatter.xsection(krill.shape, 120e3, 90, sound.speed.water, density.water)
TS <- 10 * log10(sigma.bs)

frequencies <- seq(20e3, 500e3, by=1e3)
TS.curve <- rep(0, length(frequencies))
for (i in 1:length(frequencies)) {
  TS.curve[i] <- 10 * log10(backscatter.xsection(krill.shape, frequencies[i], 90, 
                                                 sound.speed.water, density.water))
}
plot(frequencies, TS.curve, ty='l')

angles <- 180:0
TS.curve <- rep(0, length(angles))
for (i in 1:length(angles)) {
  TS.curve[i] <- 10 * log10(backscatter.xsection(krill.shape, 120, angles[i], 
                                                 sound.speed.water, density.water))
}
plot(angles, TS.curve, ty='l')
