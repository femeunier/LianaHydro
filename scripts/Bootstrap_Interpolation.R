rm(list = ls())

library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)
library(ggplot2)

FN <- c("weibull","sigmoidal","polynomial","polynomial2")


data.file <- file.path(getwd(),"data","DataOptical_Cissus_Interpolate.csv")
data <- read.csv(data.file,header = TRUE,sep=";",stringsAsFactors=FALSE,
                 na.strings=c("","NA"))

#Inteproleren, aanpassen aan aantal takken

psi <- na.omit(data[1])
PLC <- na.omit(data[2])
tak1 <- data.frame(psi,PLC,curve = 1) %>% rename(psi = WP1,PLC = Cum1)
psi <- na.omit(data[3])
PLC <- na.omit(data[4])
tak2 <- data.frame(psi,PLC,curve = 2) %>% rename(psi = WP2,PLC = Cum2)
psi <- na.omit(data[5])
PLC <- na.omit(data[6])
tak3 <- data.frame(psi,PLC,curve = 3) %>% rename(psi = WP3,PLC = Cum3)
psi <- na.omit(data[7])
PLC <- na.omit(data[8])
tak4 <- data.frame(psi,PLC,curve = 4) %>% rename(psi = WP4,PLC = Cum4)
psi <- na.omit(data[9])
PLC <- na.omit(data[10])
tak5 <- data.frame(psi,PLC,curve = 5) %>% rename(psi = WP5,PLC = Cum5)

curve.select <- seq(1,5)
data=rbind(tak1,tak2,tak3,tak4,tak5) %>% filter(curve %in% curve.select) # ###########VERANDER


library(pracma)
tak1_all <- rbind(c(0,0,3),tak1,c(min(tak1$psi)*2,100,3))
tak1_interpol <- list(y = seq(min(tak1_all$psi),0,length.out = 1000),
                      x = interp1(tak1_all$psi,tak1_all$PLC,seq(min(tak1_all$psi),0,length.out = 1000)))

tak2_all <- rbind(c(0,0,3),tak2,c(min(tak2$psi)*2,100,3))
tak2_interpol <- list(y = seq(min(tak2_all$psi),0,length.out = 1000),
                      x = interp1(tak2_all$psi,tak2_all$PLC,seq(min(tak2_all$psi),0,length.out = 1000)))

tak3_all <- rbind(c(0,0,3),tak3,c(min(tak3$psi)*2,100,3))
tak3_interpol <- list(y = seq(min(tak3_all$psi),0,length.out = 1000),
                      x = interp1(tak3_all$psi,tak3_all$PLC,seq(min(tak3_all$psi),0,length.out = 1000)))

tak4_all <- rbind(c(0,0,3),tak4,c(min(tak4$psi)*2,100,3))
tak4_interpol <- list(y = seq(min(tak4_all$psi),0,length.out = 1000),
                      x = interp1(tak4_all$psi,tak4_all$PLC,seq(min(tak4_all$psi),0,length.out = 1000)))

tak5_all <- rbind(c(0,0,3),tak5,c(min(tak3$psi)*2,100,3))
tak5_interpol <- list(y = seq(min(tak5_all$psi),0,length.out = 1000),
                      x = interp1(tak5_all$psi,tak5_all$PLC,seq(min(tak5_all$psi),0,length.out = 1000)))


# tak1_interpol=approx(tak1$PLC, tak1$psi, method="linear", n=200)
# tak2_interpol=approx(tak2$PLC, tak2$psi, method="linear", n=200)
# tak3_interpol=approx(tak3_all$PLC, tak3_all$psi, method="linear", n=200)
# tak4_interpol=approx(tak4$PLC, tak4$psi, method="linear", n=200)
# tak5_interpol=approx(tak5$PLC, tak5$psi, method="linear", n=200)

samen1<-cbind(tak1_interpol[["y"]],tak1_interpol[["x"]],1)
samen2<-cbind(tak2_interpol[["y"]],tak2_interpol[["x"]],2)
samen3<-cbind(tak3_interpol[["y"]],tak3_interpol[["x"]],3)
samen4<-cbind(tak4_interpol[["y"]],tak4_interpol[["x"]],4)
samen5<-cbind(tak5_interpol[["y"]],tak5_interpol[["x"]],5)

samen_tot<-rbind(samen1,samen2,samen3,samen4,samen5) ###########VERANDER
samen_tot<- data.frame(samen_tot) %>% rename(psi = X1,PLC = X2,curve = X3) %>% filter(curve %in% curve.select)

#Bootstrap

Nbootstrap = 250

bootstrap <- data.summary <- data.frame()
Npoints <- nrow(samen_tot)


psi_all <- seq(-9.13,0,length.out = 250) ###############meest negatieve WP

for (i in seq(1,Nbootstrap)){

  print(i/Nbootstrap)

  # Liana
  sample <- sample.int(Npoints, size = Npoints, replace = TRUE)

  liana_sample <- samen_tot[sample,]

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
  scale_x_continuous(expand = c(0,0),limits = c(-2.5,0)) +
  scale_y_continuous(limits = c(0,100),expand = c(0.01,0.01)) +
  geom_point(data = data,aes(x = psi, y = PLC)) +
  theme_bw() +
  labs(y = "PLC", x = expression(psi~'[MPa]') )

PLC<-best.modelL[["invert"]]

setwd("~/master thesis hydraulic properties lianas/Verwerking/PLC AE")
write.table(PLC,file=paste(c('CissusOptE_P50.txt'),collapse = ''),row.names=F,col.names=T)

