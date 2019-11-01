rm(list = ls())

library(LianaHydro)

lambda <- 1.5
k <- 2

psi <- seq(-0.001,-10,length.out = 1000)
PLC <- weibull(psi,lambda = lambda, k = k)
plot(psi,PLC,type='l')

psi_prim <- invert.weibull(PLC,lambda = lambda, k = k)
plot(psi_prim-psi_prim,type='l')

Slopes <-slope.weibull(c(12,50,88),lambda = lambda, k = k)

#########################################################
a <- 2
b <- -2

psi <- seq(-0.001,-10,length.out = 1000)
PLC <- sigmoidal(psi,a = a, b = b)
plot(psi,PLC,type='l')

psi_prim <- invert.sigmoidal(PLC,a = a, b = b)
plot(psi_prim-psi_prim,type='l')

Slopes <- slope.sigmoidal(c(12,50,88),a = a, b = b)

#########################################################
a <- 2
b <- -2

psi <- seq(-0.001,-10,length.out = 1000)
PLC <- polynomial(psi,a = a, b = b)
plot(psi,PLC,type='l')

psi_prim <- invert.polynomial(PLC,a = a, b = b)
plot(psi_prim-psi_prim,type='l')

Slopes <- slope.polynomial(c(12,50,88),a = a, b = b)


#########################################################
a <- 2
b <- -2

psi <- seq(-0.001,-10,length.out = 1000)
PLC <- polynomial2(psi,a = a, b = b)
plot(psi,PLC,type='l')

psi_prim <- invert.polynomial2(PLC,a = a, b = b)
plot(psi_prim-psi_prim,type='l')

Slopes <- slope.polynomial2(c(12,50,88),a = a, b = b)

