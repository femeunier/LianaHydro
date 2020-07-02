rm(list= ls())

library(LianaHydro)

fun2test <- "weibull"
default <- default.model.params(fun2test)

params <- c(default[[1]]*3,default[[2]]*2)

functionSopt <- match.fun(paste0("slope.",fun2test))
functionIopt <- match.fun(paste0("invert.",fun2test))
P50 <-  do.call(functionIopt,list(50,params[1],params[2]))
ax50 <-  do.call(functionSopt,list(50,params[1],params[2]))

params
obs2params(ax50,P50,funbest = fun2test)
