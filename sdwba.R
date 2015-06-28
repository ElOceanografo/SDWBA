library(numDeriv)

# convenience functions for vectors and trig
norm <- function(x) sqrt(sum(x^2))
dot <- function(x, y) sum(x * y)
rad2deg <- function(x) x * 180 / pi
deg2rad <- function(x) x * pi / 180



Scatterer <- function(x, y, z, a, density, sound.speed) {
  scatterer <- data.frame(x=x, y=y, z=z, a=a, density=density,
                          sound.speed=sound.speed)
  class(scatterer) <- append(class(scatterer), "Scatterer")
  return(scatterer)
}

add.density <- function(df, density.water) {
  if (! "h" %in% names(df)) {
    stop("Need a column named 'h' of density contrasts to calculate densities in animal")
  }
  if (is.null(density.water)) {
    stop("No value for density.water provided")
  }
  df$density <- df$h * density.water
  return(df)
}

add.sound.speed <- function(df, sound.speed.water) {
  if (! "g" %in% names(df)) {
    stop("Need a column named 'g' of sound speed contrasts to calculate sound speeds in animal")
  }
  if (is.null(sound.speed.water)) {
    stop("No value for sound.speed.water provided")
  }
  df$sound.speed <- df$g * sound.speed.water
  return(df)
}

load.scatterer <- function(filename, sound.speed.water=NULL,
                           density.water=NULL, ...) {
  df <- read.csv(filename, ...)
  if (! "density" %in% names(df)) {
    df <- add.density(df, density.water)
  }
  if (! "sound.speed" %in% names(df)) {
    df <- add.sound.speed(df, sound.speed.water)
  }
  return(Scatterer(x=df$x, y=df$y, z=df$z, a=df$a, density=df$density, sound.speed=df$sound.speed))
}

save.scatterer <- function(scatterer, ...)  UseMethod("save.scatterer", scatterer)

save.scatterer.Scatterer <- function(scatterer, filename, ...) {
  write.csv(scatterer, filename, ...)
}

rotate <- function(scatterer, ...) UseMethod("rotate", scatterer)

rotate.Scatterer <- function(scatterer, roll=0, tilt=0, yaw=0) {
  tilt <- deg2rad(tilt)
  roll <- deg2rad(roll)
  yaw <- deg2rad(yaw)
  Rx <- matrix(c(1, 0, 0, 0, cos(roll), sin(roll), 0, -sin(roll), cos(roll)), nrow=3)
  Ry <- matrix(c(cos(tilt), 0, -sin(tilt), 0, 1, 0, sin(tilt), 0, cos(tilt)), nrow=3)
  Rz <- matrix(c(cos(yaw), sin(yaw), 0,  -sin(yaw), cos(yaw), 0, 0, 0, 1), nrow=3)
  R <- Rz %*% Ry %*% Rx
  for (i in 1:nrow(scatterer)) {
    scatterer[i, c('x', 'y', 'z')] <- R %*% unlist(scatterer[i, c('x', 'y', 'z')])
  }
  return(scatterer)
}
  


form.function <- function(scatterer, k, sound.speed, density, phase.sd=0) {
  n <- nrow(scatterer)
  f.bs <- 0 + 0i
  gg <- scatterer$sound.speed / sound.speed
  hh <- scatterer$density / density
  for (j in 1:(n - 1)) {
    r1 <- c(scatterer$x[j], scatterer$y[j], scatterer$z[j])
    r2 <- c(scatterer$x[j + 1], scatterer$y[j + 1], scatterer$z[j + 1])
    alphatilt <- acos(dot(k, (r2 - r1)) / (norm(k) * norm(r2 - r1)))
    betatilt <- abs(alphatilt - pi/2)
    
    # Define the function to integrate over this segment
    integrand <- function(s, real=TRUE) {
      rx <- s * (r2[1] - r1[1]) + r1[1]
      ry <- s * (r2[2] - r1[2]) + r1[2]
      rz <- s * (r2[3] - r1[3]) + r1[3]
      r <- c(rx, ry, rz)
      a <- s * (scatterer$a[j + 1] - scatterer$a[j]) + scatterer$a[j]
      h <- s * (hh[j + 1] - hh[j]) + hh[j]
      g <- s * (gg[j + 1] - gg[j]) + gg[j]
      gamgam <- 1 / (g * h^2) + 1 / g - 2

      if (abs(betatilt) < 1e-10) {
        bessy <- norm(k) * a / h
      } else {
        arg <- 2 * norm(k) * a / h * cos(betatilt)
        bessy <- besselJ(arg, 1) / cos(betatilt)
      }
      result <- norm(k) / 4 * gamgam * exp(2i * dot(k, r) / h) * a * 
        bessy * norm(r2 - r1)
      if (real) {
        return(Re(result))
      } else {
        return(Im(result))
      }
    }
    integrand <- Vectorize(integrand)
    
    real.part <- integrate(integrand, 0, 1, real=TRUE)$value 
    imag.part <- integrate(integrand, 0, 1, real=FALSE)$value 
    segment.f.bs <- (real.part + imag.part) * exp(1i * rnorm(1, 0, phase.sd))
    f.bs <- f.bs + segment.f.bs
  }
  return(f.bs)
}


form.function.continuous <- function(scatterer, k, sound.speed, density,
                                     phase.sd=0) {
  # Note: phase.sd has a different meaning here than in the 
  dx <- diff(scatterer$x)
  dy <- diff(scatterer$y)
  dz <- diff(scatterer$z)
  ss <- c(0, cumsum(sqrt(dx^2 + dy^2 + dz^2)))
  g <- scatterer$sound.speed / sound.speed
  h <- scatterer$density / density
  
  # Defining continuous interpolation functions
  x_fun <- splinefun(ss, scatterer$x)
  y_fun <- splinefun(ss, scatterer$y)
  z_fun <- splinefun(ss, scatterer$z)
  a_fun <- splinefun(ss, scatterer$a)
  g_fun <- splinefun(ss, g)
  h_fun <- splinefun(ss, h)
  r_fun <- function(s) c(x_fun(s), y_fun(s), z_fun(s))
  local_tangent <- function(s) jacobian(r_fun, s)
  dot <- function(x, y) sum(x * y)
  vecnorm <- function(x) sqrt(sum(x^2))
  
  # function to integrate along the length of the animal
  integrand <- function(s, k) {
    loc_tan <- local_tangent(s)
    beta <- acos(dot(k, loc_tan) / (vecnorm(k) * vecnorm(loc_tan)))
    beta <- abs(beta - pi/2)
    (g_fun(s) - h_fun(s)) * exp(2i * dot(k, r_fun(s))) * 
      besselJ(2 * vecnorm(k) * a_fun(s) * cos(beta), 1) / cos(beta)
  }
  integrand.real <- Vectorize(function(s, k) Re(integrand(s, k)), vectorize.args=c("s"))
  integrand.imag <- Vectorize(function(s, k) Im(integrand(s, k)), vectorize.args=c("s"))
  
  f.bs.real <- integrate(integrand.real, 0, max(ss), k=k)$value
  f.bs.imag <- integrate(integrand.imag, 0, max(ss), k=k)$value
  f.bs <- (f.bs.real + f.bs.imag) * exp(1i * rnorm(1, 0, phase.sd))
  return(f.bs)
}


backscatter.xsection <- function(scatterer, k, sound.speed, density, phase.sd=0,
                                 method=c("discrete", "continuous")) {
  method <- match.arg(method)
  if (method == "discrete") {
    f.bs <- form.function(scatterer, k, sound.speed, density, phase.sd)
  } else if (method == "continuous") {
    f.bs <- form.function.continuous(scatterer, k, sound.speed, density, phase.sd)
  } else {
    stop("method must be either 'discrete' or 'continuous'")
  }
  return(abs(f.bs)^2)
}


TS <- function(scatterer, k, sound.speed, density, phase.sd=0,
               method=c("discrete", "continuous")) {
  sigma.bs <- backscatter.xsection(scatterer, k, sound.speed, density, phase.sd, method)
  return(10 * log10(sigma.bs))
}

frequency.spectrum <- function(scatterer, freq.start, freq.stop, sound.speed, 
                               density, nfreq=100, tilt.angle=0, ...) {
  freqs <- seq(freq.start, freq.stop, length.out=nfreq)
  sigma.bs <- rep(0, length(freqs))
  for (i in 1:nfreq) {
    k <- c(0, 0, 2 * pi * freqs[i] / sound.speed)
    sigma.bs[i] <- backscatter.xsection(scatterer, k, sound.speed, density, ...)
  }
  return(list(freq=freqs, sigma.bs=sigma.bs, TS=10*log10(sigma.bs)))
}


tilt.spectrum <- function(scatterer, angle.start, angle.stop, freq, sound.speed,
                          density, nangle=100, ...) {
  angles <- seq(angle.start, angle.stop, length.out=nangle) # degrees
  sigma.bs <- rep(0, length(angles))
  for (i in 1:length(angles)) {
    k <- c(0, 0, 2 * pi * freq / sound.speed)
    sigma.bs[i] <- backscatter.xsection(rotate(scatterer, tilt=angles[i]), k, 
                                        sound.speed, density, ...)
  }
  return(list(angle=angles, sigma.bs=sigma.bs, TS=10 * log10(sigma.bs)))
}
