rm(list = ls())

library(LianaHydro)

lambda <- 1.5
k <- 2

psi <- seq(-0.001,-10,length.out = 1000)
PLC <- weibull(psi,param = c(lambda,k))
plot(psi,PLC,type='l')

psi_prim <- invert.weibull(PLC,param=c(lambda,k))
plot(psi_prim-psi_prim,type='l')

Slopes <-slope.weibull(c(12,50,88),param=c(lambda,k))
