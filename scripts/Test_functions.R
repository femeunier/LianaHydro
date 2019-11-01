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
PLC <- sigmoidal(psi,param = c(a,b))
plot(psi,PLC,type='l')

psi_prim <- invert.sigmoidal(PLC,param = c(a,b))
plot(psi_prim-psi_prim,type='l')

Slopes <- slope.sigmoidal(c(12,50,88),param = c(a,b))
