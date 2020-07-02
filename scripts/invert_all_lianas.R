rm(list = ls())

library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)
library(ggplot2)
library(reshape2)
library(zoo)

file <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
dataVC <- read.csv(file)

FN <- c("weibull","sigmoidal","polynomial","polynomial2","cumnorm")
Ps <- c(12,50,88)

Cols <- c(rgb(0,0,139/255),rgb(0.10,0.50,0.00))
data.files <- c("./data/treerawdata.csv","./data/lianarawdata.csv")
GF <- c("Tree","Liana")

invert.PLC <- summary <- data.frame()

data.all <- data.frame()

psi <- seq(-10,0,length.out = 100)
maxplant <- 0
lianas <- trees <- c()
for (ifile in seq(data.files)){

  data <- read.csv(data.files[ifile],header = TRUE) %>% mutate(psi = - abs(psi))

  data.all <- rbind(data.all,
                    data %>% dplyr::select(Id,psi,PLC) %>% mutate(GF = GF[ifile]))

  plants <- sort(unique(data$Id))


  best.models <- list()
  for (plant in plants){
    data.plant <- data %>% filter(Id == plant)
    models <- opt.data(data = data.plant,
                       function.names = FN)
    models <- add.properties(models,x = Ps)
    best.model <- find.best.model(models)[[1]]
    best.models[[as.character(plant)]] <- best.model

    summary <- rbind(summary,
                     data.frame(Id = plant+maxplant,
                                GF = GF[ifile],
                                model = best.model$name,
                                P50 = best.model$invert[which(names(best.model$invert) == "P50")],
                                ax50 = best.model$slopes[which(names(best.model$slopes) == "ax50")],
                                RMSE = best.model$RMSE,
                                r2 = best.model$r.squared))

    functionopt <- match.fun(best.model$name)
    PLC <- do.call(functionopt,list(psi,best.model$best.params[1],best.model$best.params[2]))



    if (GF[ifile] == "Tree"){
      trees <- rbind(trees,PLC)
      k = mean(dataVC %>% filter(GrowthForm == "Tree" & Id == plant) %>% pull(ksat),na.rm=TRUE)
    } else {
      lianas <- rbind(lianas,PLC)
      k = mean(dataVC %>% filter(GrowthForm == "Tree" & Id == plant) %>% pull(ksat),na.rm=TRUE)
    }

    invert.PLC <- rbind(invert.PLC,
                        data.frame(psi = psi,
                                   PLC = PLC,
                                   Id = maxplant+plant,
                                   GF = GF[ifile],
                                   k = k*(1- PLC/100)))

  }
  maxplant <- max(invert.PLC %>% pull(Id))
  print(maxplant)
}

ggplot(data = data.all,
       aes(x = psi, y = PLC, color = as.factor(GF),group = Id)) +
  geom_point() +
  geom_line(data = invert.PLC,aes()) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  theme_bw()


data.sum <- invert.PLC %>% group_by(GF,psi) %>% summarise(PLC_m = mean(PLC,na.rm = TRUE),
                                                 PLC_low = quantile(PLC,0.025,na.rm = TRUE),
                                                 PLC_high = quantile(PLC,0.975,na.rm = TRUE),
                                                 k_m = mean(k,na.rm = TRUE),
                                                 k_low = quantile(k,0.025,na.rm = TRUE),
                                                 k_high = quantile(k,0.975,na.rm = TRUE))


# bootstraping

Nbootstrap = 250
Nliana <- nrow(lianas)
Ntree <- nrow(trees)

ksat_liana <- dataVC %>% filter(!is.na(ksat) & GrowthForm == "Liana")%>%pull(ksat)
Nliana_K <- length(ksat_liana)

ksat_tree <- dataVC %>% filter(!is.na(ksat) & GrowthForm == "Tree")%>%pull(ksat)
Ntree_K <- length(ksat_tree)

bootstrap <- data.summary <- data.frame()

for (i in seq(1,Nbootstrap)){

  print(i/Nbootstrap)

  # Liana
  sample <- sample.int(Nliana, size = Nliana, replace = TRUE)

  liana_sample <- lianas[sample,]
  data.liana <- data.frame(psi = rep(psi,Nliana),
                           PLC = as.vector(t(liana_sample)))
  models <- opt.data(data = data.liana,
                     function.names = c("weibull","sigmoidal","polynomial","polynomial2","cumnorm"))
  models <- add.properties(models,x = 50)
  best.modelL <- find.best.model(models)[[1]]

  sampleK <- sample.int(Nliana_K, size = Nliana_K, replace = TRUE)
  KsampleL <- mean(ksat_liana[sampleK])

  PLC_L <- best.modelL$PLC.predict.all
  K_L <- KsampleL*(1-PLC_L/100)

  bootstrap <- rbind(bootstrap,
                     data.frame(id = i,
                                psi = best.modelL$psi.all,
                                PLC = PLC_L,
                                k =  K_L,GF = "Liana"))

  # Tree
  sample <- sample.int(Ntree, size = Ntree, replace = TRUE)

  tree_sample <- trees[sample,]
  data.tree  <- data.frame(psi = rep(psi,Ntree),
                           PLC = as.vector(t(tree_sample)))
  models <- opt.data(data = data.tree,
                     function.names = c("weibull","sigmoidal","polynomial","polynomial2"))
  models <- add.properties(models,x = 50)
  best.modelT <- find.best.model(models)[[1]]

  sampleK <- sample.int(Ntree_K, size = Ntree_K, replace = TRUE)
  KsampleT <- mean(ksat_tree[sampleK])

  PLC_T <- best.modelT$PLC.predict.all
  K_T <- KsampleT*(1-PLC_T/100)

  bootstrap <- rbind(bootstrap,
                     data.frame(id = i,
                                psi = best.modelT$psi.all,
                                PLC = PLC_T,
                                k =  K_T,
                                GF = "Tree"))

  # Summary
  data.summary <- rbind(data.summary,
                        data.frame(P50 = best.modelL$invert,
                                   ax50 = best.modelL$slopes,
                                   GF = "Liana"),
                        data.frame(P50 = best.modelT$invert,
                                   ax50 = best.modelT$slopes,
                                   GF = "Tree"))
}


N = 50
signif.all <- data.frame()
for (i in seq(1,N)){
  print(i)
  signif <- bootstrap %>% group_by(psi,GF) %>% sample_n(5) %>% ungroup() %>% group_by(psi) %>% summarise(PLC = kruskal.test(formula = PLC ~ GF)$p.value,
                                                                                                          k  = kruskal.test(formula = k ~ GF)$p.value) %>% mutate(num = i)
  signif.all <- rbind(signif.all,
                      signif)
}

signif <- signif.all %>% group_by(psi) %>% summarise(alpha_PLC = mean(PLC),
                                                     alpha_k = mean(k)) %>% mutate(alpha_PLC = rollapply(alpha_PLC,10,mean,na.rm=TRUE,partial=TRUE),
                                                                                   alpha_k = rollapply(alpha_k,10,mean,na.rm=TRUE,partial=TRUE))

bootstrap_sum <- bootstrap %>% filter(k > 0)  %>% group_by(GF,psi) %>% summarise(PLC_m = mean(PLC),
                                                                                 PLC_low = quantile(PLC,0.025),
                                                                                 PLC_high = quantile(PLC,0.975),
                                                                                 k_m = mean(k),
                                                                                 k_low = quantile(k,0.025),
                                                                                 k_high = quantile(k,0.975))

pos <- bootstrap_sum %>% group_by(GF) %>% summarise(P50low = psi[which.min(abs(PLC_low - 50))],
                                                        P50high = psi[which.min(abs(PLC_high - 50))],
                                                        P50m = psi[which.min(abs(PLC_m - 50))])


P50l <- pos %>% filter(GF == "Liana")
P50t <- pos %>% filter(GF == "Tree")

slopes <- data.summary %>% group_by(GF) %>% summarise(ax50m = mean(ax50))
intercept <- c(50-slopes$ax50m[1]*P50l$P50m,50-slopes$ax50m[2]*P50t$P50m)
psi_x <- c(0.5,1.5)
y = slopes$ax50m[1]*psi_x*c(P50l$P50m)+intercept[1]

psi_x2 <- c(0.6,1.4)
y2 = slopes$ax50m[2]*psi_x2*c(P50t$P50m)+intercept[2]

w = 1; Top = 100 ; bottom = 1
# PLC curves
ggplot() +
  geom_ribbon(data = bootstrap_sum,aes(x = psi,color = as.factor(GF),
                                       fill = as.factor(GF),ymin = PLC_low,ymax = PLC_high),alpha = 0.1,colour = NA) +
  geom_ribbon(data = data.frame(x = c(P50l$P50low,P50l$P50high),
                                ymin = bottom -c(w,w),
                                ymax = bottom +c(w,w)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[1], colour = NA,alpha = 0.1) +
  geom_ribbon(data = data.frame(x = c(P50l$P50m,P50l$P50m)*c(0.99,1.01),
                                ymin = bottom -c(w,w),
                                ymax = bottom +c(w,w)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[1], colour = NA,alpha = 1) +
  geom_ribbon(data = data.frame(x = c(P50t$P50low,P50t$P50high),
                                ymin = bottom -c(w,w),
                                ymax = bottom +c(w,w)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[2], colour = NA,alpha = 0.1) +
  geom_ribbon(data = data.frame(x = c(P50t$P50m,P50t$P50m)*c(0.995,1.01),
                                ymin = bottom -c(w,w),
                                ymax = bottom +c(w,w)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[2], colour = NA,alpha = 1) +
  geom_ribbon(data = data.frame(x = signif$psi[signif$alpha_PLC<0.01],
                                ymin = Top -w,
                                ymax = Top + w),aes(x = x,ymin = ymin, ymax = ymax), fill = "darkgrey", colour = NA,alpha = 0.5) +
  geom_ribbon(data = data.frame(x = signif$psi[signif$alpha_PLC<0.05],
                                ymin = Top -w,
                                ymax = Top + w),aes(x = x,ymin = ymin, ymax = ymax), fill = "lightgrey", colour = NA,alpha = 0.5) +
  geom_segment(aes(x = psi_x[1]*P50l$P50m,xend = psi_x[2]*P50l$P50m,
                   y = y[1], yend = y[2]), colour = Cols[1],linetype = 2) +
  geom_segment(aes(x = psi_x2[1]*P50t$P50m,xend = psi_x2[2]*P50t$P50m,
                   y = y2[1], yend = y2[2]), colour = Cols[2],linetype = 2) +
 geom_line(data = bootstrap_sum,aes(x = psi,y = PLC_m,color = as.factor(GF))) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,102),expand = c(0.0,0.0)) +
  theme_bw() + theme(legend.position = "none")

ggplot(data = dataVC, aes(x = GrowthForm,y = ksat,
                              fill = as.factor(GrowthForm)))+
  geom_boxplot(alpha = 0.3) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  scale_y_log10() +
  labs(x = "") +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_blank())

PmdL <- quantile(dataVC %>% filter(!is.na(Pmd) & GrowthForm == "Liana") %>% pull(Pmd),c(0.025,0.5,0.975))
PmdT <- quantile(dataVC %>% filter(!is.na(Pmd) & GrowthForm == "Tree") %>% pull(Pmd),c(0.025,0.5,0.975))


# k curves
Top = 35 ; w = 4 ; bottom  = 0.003 ; wbot = 0.0002 ; bottom2 = bottom  - 2*wbot; wbot2 = 0.00001
ggplot() +
  geom_ribbon(data = bootstrap_sum,aes(x = psi,color = as.factor(GF),
                                       fill = as.factor(GF),ymin = k_low,ymax = k_high),alpha = 0.1,colour = NA) +
  geom_line(data = bootstrap_sum,aes(x = psi,y = k_m,color = as.factor(GF))) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  geom_ribbon(data = data.frame(x = signif$psi[signif$alpha_k<0.01],
                                ymin = Top -w,
                                ymax = Top + w),aes(x = x,ymin = ymin, ymax = ymax), fill = "darkgrey", colour = NA,alpha = 0.5) +
  geom_ribbon(data = data.frame(x = signif$psi[signif$alpha_k<0.05],
                                ymin = Top -w,
                                ymax = Top + w),aes(x = x,ymin = ymin, ymax = ymax), fill = "lightgrey", colour = NA,alpha = 0.5) +
  geom_ribbon(data = data.frame(x = c(PmdL[1],PmdL[3]),
                                ymin = bottom -c(wbot,wbot),
                                ymax = bottom +c(wbot,wbot)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[1], colour = NA,alpha = 0.1) +
  geom_ribbon(data = data.frame(x = c(PmdL[2],PmdL[2])*c(0.995,1.01),
                                ymin = bottom -c(wbot,wbot),
                                ymax = bottom +c(wbot,wbot)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[1], colour = NA,alpha = 1) +
  geom_ribbon(data = data.frame(x = c(PmdT[1],PmdT[3]),
                                ymin = bottom2 -c(wbot,wbot),
                                ymax = bottom2 +c(wbot,wbot)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[2], colour = NA,alpha = 0.1) +
  geom_ribbon(data = data.frame(x = c(PmdT[2],PmdT[2])*c(0.99,1.01),
                                ymin = bottom2 -c(wbot,wbot),
                                ymax = bottom2 +c(wbot,wbot)),aes(x = x,ymin = ymin, ymax = ymax), fill = Cols[2], colour = NA,alpha = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0.01,0.01)) +
  theme_bw() + theme(legend.position = "none")



