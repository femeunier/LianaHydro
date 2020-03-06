rm(list = ls())

library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)
library(ggplot2)

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

# # PLC curves no bootstraping
# ggplot(data = data.sum,aes(x = psi,y = PLC_m,color = as.factor(GF),
#                                 fill = as.factor(GF),ymin = PLC_low,ymax = PLC_high)) +
#   geom_ribbon(alpha = 0.1,colour = NA) +
#   geom_line() +
#   scale_color_manual(values = rev(Cols)) +
#   scale_fill_manual(values = rev(Cols)) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(limits = c(0,100),expand = c(0.01,0.01)) +
#   theme_bw()
#
# ggplot(data = data.sum,aes(x = psi,y = k_m,color = as.factor(GF),
#                            fill = as.factor(GF),ymin = k_low,ymax = k_high)) +
#   geom_ribbon(alpha = 0.1,colour = NA) +
#   geom_line() +
#   scale_color_manual(values = rev(Cols)) +
#   scale_fill_manual(values = rev(Cols)) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_log10(expand = c(0.01,0.01)) +
#   theme_bw()



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
  signif <- bootstrap %>% group_by(psi,GF) %>% sample_n(15) %>% ungroup() %>% group_by(psi) %>% summarise(PLC = kruskal.test(formula = PLC ~ GF)$p.value,
                                                                                                          k  = kruskal.test(formula = k ~ GF)$p.value) %>% mutate(num = i)
  signif.all <- rbind(signif.all,
                      signif)
}

signif <- signif.all %>% group_by(psi) %>% summarise(alpha_PLC = mean(PLC),
                                                     alpha_k = mean(k))

bootstrap_sum <- bootstrap %>% filter(k > 0)  %>% group_by(GF,psi) %>% summarise(PLC_m = mean(PLC),
                                                                                 PLC_low = quantile(PLC,0.025),
                                                                                 PLC_high = quantile(PLC,0.975),
                                                                                 k_m = mean(k),
                                                                                 k_low = quantile(k,0.025),
                                                                                 k_high = quantile(k,0.975))


# PLC curves
ggplot(data = bootstrap_sum,aes(x = psi,y = PLC_m,color = as.factor(GF),
                                fill = as.factor(GF),ymin = PLC_low,ymax = PLC_high)) +
  geom_ribbon(alpha = 0.1,colour = NA) +
  geom_line() +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,100),expand = c(0.01,0.01)) +
  theme_bw()

# k curves
ggplot(data = bootstrap_sum,aes(x = psi,y = k_m,color = as.factor(GF),
                                fill = as.factor(GF),ymin = k_low,ymax = k_high)) +
  geom_ribbon(alpha = 0.1,colour = NA) +
  geom_line() +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0.01,0.01)) +
  theme_bw()


quantile(dataVC %>% filter(!is.na(p50) & GrowthForm == "Liana") %>% pull(p50),c(0.025,0.5,0.975))
quantile(dataVC %>% filter(!is.na(p50) & GrowthForm == "Tree") %>% pull(p50),c(0.025,0.5,0.975))



quantile(dataVC %>% filter(!is.na(Pmd) & GrowthForm == "Liana") %>% pull(Pmd),c(0.025,0.5,0.975))
quantile(dataVC %>% filter(!is.na(Pmd) & GrowthForm == "Tree") %>% pull(Pmd),c(0.025,0.5,0.975))
