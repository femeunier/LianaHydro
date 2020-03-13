rm(list = ls())

library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)

data.file <- "./data/lianarawdata.csv"
data <- read.csv(data.file,header = TRUE) %>% rename(psi = psi,
                                                     PLC = PLC)

data.test <- data %>% filter(Id == 5) %>% dplyr::select(psi,PLC) %>% mutate(psi = pmin(0,-abs(psi)))

All.models <-
  opt.data(data = data.test)
length(All.models)
All.models <- add.properties(All.models,x = c(12,50,88))


best.model <-
  find.best.model(All.models)

# plot.models(All.models,add=TRUE)
plot(data.test[["psi"]],data.test[["PLC"]],pch=19,xlab = "Psi",ylab = "PLC",ylim = c(0,100))
plot.models(All.models,col=c('black'),add=TRUE,highlight=TRUE)


