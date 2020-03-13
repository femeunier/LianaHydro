rm(list = ls())

library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)
library(ggplot2)

FN <- c("weibull","sigmoidal","polynomial","polynomial2")

data.file <- "./data/lianarawdata.csv"
data <- read.csv(data.file,header = TRUE) %>% rename(psi = psi,
                                                     PLC = PLC)

data.test <- data %>% filter(Id == 17) %>% dplyr::select(psi,PLC) %>% mutate(psi = pmin(0,-abs(psi)))

Nbootstrap = 250

bootstrap <- data.summary <- data.frame()
Npoints <- nrow(data.test)


psi_all <- seq(-10,0,length.out = 250)

for (i in seq(1,Nbootstrap)){

  print(i/Nbootstrap)

  # Liana
  sample <- sample.int(Npoints, size = Npoints, replace = TRUE)

  liana_sample <- data.test[sample,]

  models <- opt.data(data = liana_sample,
                     function.names = FN)
  models <- add.properties(models,x = 50)
  best.modelL <- find.best.model(models)[[1]]

  functionopt <- match.fun(best.modelL$name)
  PLC <- do.call(functionopt,list(psi_all,best.modelL$best.params[1],best.modelL$best.params[2]))

  bootstrap <- rbind(bootstrap,
                     data.frame(id = i,
                                psi = psi_all,
                                PLC = PLC))
  # Summary
  data.summary <- rbind(data.summary,
                        data.frame(P50 = best.modelL$invert,
                                   ax50 = best.modelL$slopes))
}

bootstrap_sum <- bootstrap %>% group_by(psi) %>% summarise(PLC_m = mean(PLC),
                                                           PLC_low = quantile(PLC,0.025),
                                                           PLC_high = quantile(PLC,0.975))

#P50 distribution
boxplot(data.summary$P50)

# Curve envelope
ggplot() +
  geom_ribbon(data = bootstrap_sum,aes(x = psi,ymin = PLC_low,ymax = PLC_high),alpha = 0.1,colour = NA) +
  geom_line(data = bootstrap_sum,aes(x = psi,y = PLC_m)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,100),expand = c(0.01,0.01)) +
  geom_point(data = data.test,aes(x = psi, y = PLC)) +
  theme_bw()


