rm(list = ls())

library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)

data.file <- "./data/rawdata.csv"
data <- read.csv(data.file,header = TRUE) %>% rename(psi = Psi..MPa.,
                                                     PLC = PLC....)

data.test <- data %>% filter(Liana.id == 1)

m0 <- nlsLM(data = data.test,
            PLC ~ weibull(psi,lambda,k),
            start=list(lambda = 2, k = 2), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                                 printEval = TRUE, warnOnly = TRUE))

m1 <- nlsLM(data = data.test,
            PLC ~ sigmoidal(psi,a,b),
            start=list(a = 2, b = -2), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                 printEval = TRUE, warnOnly = TRUE))

m2 <- nlsLM(data = data.test,
            PLC ~ polynomial(psi,a,b),
            start=list(a = 2, b = -2), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                             printEval = TRUE, warnOnly = TRUE))

m3 <- nlsLM(data = data.test,
            PLC ~ polynomial2(psi,a,b),
            start=list(a = 2, b = -2), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                             printEval = TRUE, warnOnly = TRUE))


best.params0 <- coef(m0)
best.params1 <- coef(m1)
best.params2 <- coef(m2)
best.params3 <- coef(m3)

RMSE <- rep(NA,4)
RMSE[1] <- weibull.comp(data.test,k = best.params0["k"],lambda = best.params0["lambda"])$RMSE
RMSE[2] <- sigmoidal.comp(data.test,a = best.params1["a"],b = best.params1["b"])$RMSE
RMSE[3] <- polynomial.comp(data.test,a = best.params2["a"],b = best.params2["b"])$RMSE
RMSE[4] <- polynomial2.comp(data.test,a = best.params3["a"],b = best.params3["b"])$RMSE

psi_extr <- extremum(data.test[["psi"]])
psi <- seq(psi_extr[1],psi_extr[2],length.out = 1000)

plot(data.test[["psi"]],data.test[["PLC"]],pch=19,xlab = "Psi",ylab = "PLC",ylim = c(0,100))
lines(psi,weibull(psi,k = best.params0["k"],lambda = best.params0["lambda"]),lty = 1)
lines(psi,sigmoidal(psi,a = best.params1["a"],b = best.params1["b"]),lty = 2)
lines(psi,polynomial(psi,a = best.params2["a"],b = best.params2["b"]),lty = 3)
lines(psi,polynomial2(psi,a = best.params3["a"],b = best.params3["b"]),lty = 4)

