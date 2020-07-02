rm(list = ls())
setwd("~/Documents/projects/LianaHydro/")
library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)
library(ggplot2)

FN <- c("weibull","sigmoidal","polynomial","polynomial2")


data.file <- file.path(getwd(),"data","dataOptical_Cissus_all2.csv")
data <- read.csv(data.file,header = TRUE,sep=";",stringsAsFactors=FALSE,
                 na.strings=c("","NA"))

psi <- na.omit(data[1])
PLC <- na.omit(data[2])
data <- data.frame(psi,PLC) %>% rename(psi = WP1,PLC = Cum1)


Nbootstrap = 250

bootstrap <- data.summary <- data.frame()
Npoints <- nrow(data)


psi_all <- seq(-1.65,0,length.out = 250)

for (i in seq(1,Nbootstrap)){

  print(i/Nbootstrap)

  # Liana
  sample <- sample.int(Npoints, size = Npoints, replace = TRUE)

  liana_sample <- data[sample,]

  models <- opt.data(data = liana_sample,
                     function.names = FN)
  models <- add.properties(models,x = c(12,50,88))
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
  geom_point(data = data,aes(x = psi, y = PLC)) +
  theme_bw() +
  labs(y = "PLC", x = "Waterpotentiaal (MPa)" )

PLC<-best.modelL[["invert"]]

setwd("~/master thesis hydraulic properties lianas/Verwerking/PLC AE")
write.table(PLC,file=paste(c('CissusOptE_P50.txt'),collapse = ''),row.names=F,col.names=T)

