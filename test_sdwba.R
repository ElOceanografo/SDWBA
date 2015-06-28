source("sdwba.R")

sound.speed.water <- 1480
density.water <- 1027

krill <- load.scatterer("data/generic_krill_McGeehee1998.csv", sound.speed.water, density.water)
krill <- krill[order(krill$x), ]
# convert form meters to mm
krill[c("x", "y", "z", "a")] <- krill[c("x", "y", "z", "a")] * 1e-3

k.mag <- 120e3 * 2 * pi / sound.speed.water
k <- c(0, 0, k.mag)
TS(krill, k, sound.speed.water, density.water)
TS(krill, k, sound.speed.water, density.water, phase.sd=pi/10)
TS(krill, k, sound.speed.water, density.water, method = "continuous")


fs <- frequency.spectrum(krill, 12e3, 400e3, sound.speed.water,
                         density.water, nfreq=120)
fs.cont <- frequency.spectrum(krill, 12e3, 400e3, sound.speed.water,
                         density.water, nfreq=120, method="continuous")
plot(TS ~ freq, fs, ty='l')
lines(TS ~ freq, fs.cont, col='blue')

krill.r <- rotate(krill, tilt=3)
plot(z ~ x, krill, ty='b')
points(z ~ x, krill.r, ty='b')

ts <- tilt.spectrum(krill, -90, 90, freq=200e3, sound.speed=sound.speed.water,
              density=density.water, nangle=180)
ts.cont <- tilt.spectrum(krill, -90, 90, freq=200e3, sound.speed=sound.speed.water,
              density=density.water, nangle=180, method="continuous")

plot(TS ~ angle, ts, ty='l')
lines(TS ~ angle, ts.cont, col='blue')
