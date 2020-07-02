rm(list = ls())

library(dplyr)
library(BayesianTools)
library(minpack.lm)
library(LianaHydro)
library(ggplot2)
library(reshape2)
library(zoo)
library(stringr)
library(ggstance)

filePV <- "/home/femeunier/Documents/projects/LianaHydro/data/PV.all.csv"
dataPV <- read.csv(filePV,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% filter(Organ == "Leaf") %>% dplyr::select(Species,GrowthForm,p0,tlp,rwc.tlp,epsil,Cft,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Organ) %>%
  mutate(kl=NA,ksat=NA,Al.As=NA,ax=NA,p50=NA,Pmd=NA,Ppd=NA,Id=NA,Function=NA)

fileVC <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
dataVC <- read.csv(fileVC,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% dplyr::select(Species,GrowthForm,kl,ksat,Al.As,ax,p50,Pmd,Ppd,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Id,Function) %>% mutate(
  p0=NA,tlp=NA,rwc.tlp=NA,epsil=NA,Cft=NA,Organ=NA)

data.all <- rbind(dataPV,dataVC)
data.all <- correct.obs(data.all)
data.all <- add.properties.data(data.all)

Cols <- c(rgb(0,0,139/255),rgb(0.10,0.50,0.00))

mLianas <- nlsLM(data = dataVC %>% filter(GrowthForm == "Liana"),
                 ax ~ a*(-p50)**b,
                 start=list(a=54.4, b=-1.17), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                    printEval = TRUE, warnOnly = TRUE))
mTrees <- nlsLM(data = dataVC %>% filter(GrowthForm == "Tree"),
                ax ~ a*(-p50)**b,
                start=list(a=54.4, b=-1.17), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                   printEval = TRUE, warnOnly = TRUE))

data_added <- dataVC %>% mutate(ax.mod = case_when(
  GrowthForm == "Liana" & is.na(ax) ~  coef(mLianas)[1]*(-p50)** coef(mLianas)[2],
  GrowthForm == "Tree" & is.na(ax) ~  coef(mTrees)[1]*(-p50)** coef(mTrees)[2],
  TRUE ~ ax))


N <- nrow(dataVC)

FN <- c("weibull","sigmoidal","polynomial","polynomial2","cumnorm")
Ps <- c(12,50,88)

psi <- seq(-10,0,length.out = 100)
lianas <- trees <- summary <- c()

for (i in seq(1,N)){

  currentId <- data_added$Id[i]
  currentGF <- data_added$GrowthForm[i]
  currentP50 <- data_added$p50[i]
  currentax <- data_added$ax[i]
  currentfun <- data_added$Function[i]

  if(is.finite(currentId)){
    file <- paste0("./data/",tolower(currentGF),"rawdata.csv")
    data <- read.csv(file,header = TRUE) %>% mutate(psi = - abs(psi)) %>% dplyr::select(Id,psi,PLC) %>% filter(Id == currentId)

    models <- opt.data(data = data, function.names = FN)
    models <- add.properties(models,x = Ps)
    best.model <- find.best.model(models)[[1]]

    current <-  data.frame(Id = i,
                           GF = currentGF,
                           model = best.model$name,
                           P50 = best.model$invert[which(names(best.model$invert) == "P50")],
                           ax50 = best.model$slopes[which(names(best.model$slopes) == "ax50")],
                           RMSE = best.model$RMSE,
                           r2 = best.model$r.squared)

    summary <- rbind(summary,
                    current)

    function2call <- match.fun(best.model$name)
    PLC <- do.call(function2call,list(psi,best.model$best.params[1],best.model$best.params[2]))

    if (currentGF == "Liana"){
      lianas <- rbind(lianas,PLC)
    } else {
      trees <- rbind(trees,PLC)
    }


  } else if (!is.na(currentP50) & !is.na(currentax) & !is.na(currentfun) & str_length(currentfun)>0){
    params <- obs2params(ax50=-abs(currentax),P50=currentP50,funbest = currentfun)

    current <-
      data.frame(Id = i,
                 GF = currentGF,
                 model = currentfun,
                 P50 = currentP50,
                 ax50 = -abs(currentax),
                 RMSE = NA,
                 r2 = NA)

    summary <- rbind(summary,current)
    function2call <- match.fun(currentfun)
    PLC <- do.call(function2call,list(psi,params[1],params[2]))

    if (currentGF == "Liana"){
      lianas <- rbind(lianas,PLC)
    } else {
      trees <- rbind(trees,PLC)
    }
  }


}

########################################################################################################################
# Bootstrap


Nbootstrap = 30
Nliana <- nrow(lianas)
Ntree <- nrow(trees)

ksat_liana <- dataVC %>% filter(!is.na(kl) & GrowthForm == "Liana")%>%pull(kl)
Nliana_K <- length(ksat_liana)

ksat_tree <- dataVC %>% filter(!is.na(kl) & GrowthForm == "Tree")%>%pull(kl)
Ntree_K <- length(ksat_tree)

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


N = 1
signif.all <- data.frame()
for (i in seq(1,N)){
  print(i)
  signif <- bootstrap %>% group_by(psi,GF) %>% ungroup() %>% group_by(psi) %>% summarise(PLC = kruskal.test(formula = PLC ~ GF)$p.value,
                                                                                                         k  = kruskal.test(formula = k ~ GF)$p.value) %>% mutate(num = i)
  signif.all <- rbind(signif.all,
                      signif)
}

signif <- signif.all %>% group_by(psi) %>% summarise(alpha_PLC = mean(PLC),
                                                     alpha_k = mean(k)) %>% mutate(alpha_PLC = rollapply(alpha_PLC,100,mean,na.rm=TRUE,partial=TRUE),
                                                                                   alpha_k = rollapply(alpha_k,100,mean,na.rm=TRUE,partial=TRUE))

bootstrap_sum <- bootstrap %>% filter(k > 0)  %>% group_by(GF,psi) %>% summarise(PLC_m = mean(PLC),
                                                                                 PLC_low = quantile(PLC,0.025),
                                                                                 PLC_high = quantile(PLC,0.975),
                                                                                 k_m = mean(k),
                                                                                 k_low = quantile(k,0.025),
                                                                                 k_high = quantile(k,0.975))

pos <- bootstrap_sum %>% group_by(GF) %>% summarise(P50low = psi[which.min(abs(PLC_low - 50))],
                                                    P50high = psi[which.min(abs(PLC_high - 50))],
                                                    P50m = psi[which.min(abs(PLC_m - 50))])

P50bootstrap <- bootstrap %>% group_by(GF,id) %>% filter((abs(PLC - 50)) == min(abs(PLC - 50)))

P50l <- pos %>% filter(GF == "Liana")
P50t <- pos %>% filter(GF == "Tree")

slopes <- data.summary %>% group_by(GF) %>% summarise(ax50m = mean(ax50))
intercept <- c(50-slopes$ax50m[1]*P50l$P50m,50-slopes$ax50m[2]*P50t$P50m)
psi_x <- c(0.5,1.5)
y = slopes$ax50m[1]*psi_x*c(P50l$P50m)+intercept[1]

psi_x2 <- c(0.6,1.4)
y2 = slopes$ax50m[2]*psi_x2*c(P50t$P50m)+intercept[2]

w = 0.8; Top = 101.1 ; bottom = 1
# PLC curves
ggplot() +
  geom_ribbon(data = bootstrap_sum,aes(x = psi,color = as.factor(GF),
                                       fill = as.factor(GF),ymin = PLC_low,ymax = PLC_high),alpha = 0.1,colour = NA) +
  geom_boxploth(data = P50bootstrap %>% filter(GF == "Liana"),aes(x = psi,y= 3,fill = GF),alpha = 0.2,width = 5,outlier.shape = NA) +
  geom_boxploth(data = P50bootstrap %>% filter(GF == "Tree"),aes(x = psi,y= 3,fill = GF),alpha = 0.2,width = 5,outlier.shape = NA) +
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
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 16)) + labs(x = "",y = "")

ggsave(plot = last_plot(),
       filename = "./Figures/VC.png", dpi = 300, width = 15, height = 8)

# # PLC curves
# ggplot() +
#   geom_ribbon(data = bootstrap_sum,aes(x = psi,color = as.factor(GF),
#                                        fill = as.factor(GF),ymin = PLC_low,ymax = PLC_high),alpha = 0.1,colour = NA) +
#   geom_boxploth(data = data.all %>% filter(GrowthForm == "Liana"),aes(x = p50,y= 6,fill = GrowthForm),alpha = 0.2,width = 3,outlier.shape = NA) +
#   geom_boxploth(data = data.all %>% filter(GrowthForm == "Tree"),aes(x = p50,y= 3,fill = GrowthForm),alpha = 0.2,width = 3,outlier.shape = NA) +
#   geom_ribbon(data = data.frame(x = signif$psi[signif$alpha_PLC<0.01],
#                                 ymin = Top -w,
#                                 ymax = Top + w),aes(x = x,ymin = ymin, ymax = ymax), fill = "darkgrey", colour = NA,alpha = 0.5) +
#   geom_ribbon(data = data.frame(x = signif$psi[signif$alpha_PLC<0.05],
#                                 ymin = Top -w,
#                                 ymax = Top + w),aes(x = x,ymin = ymin, ymax = ymax), fill = "lightgrey", colour = NA,alpha = 0.5) +
#   geom_segment(aes(x = psi_x[1]*P50l$P50m,xend = psi_x[2]*P50l$P50m,
#                    y = y[1], yend = y[2]), colour = Cols[1],linetype = 2) +
#   geom_segment(aes(x = psi_x2[1]*P50t$P50m,xend = psi_x2[2]*P50t$P50m,
#                    y = y2[1], yend = y2[2]), colour = Cols[2],linetype = 2) +
#   geom_line(data = bootstrap_sum,aes(x = psi,y = PLC_m,color = as.factor(GF))) +
#   scale_color_manual(values = Cols) +
#   scale_fill_manual(values = Cols) +
#   scale_x_continuous(expand = c(0,0)) +
#   scale_y_continuous(limits = c(0,102),expand = c(0.0,0.0)) +
#   theme_bw() + theme(legend.position = "none") + labs(x = "",y = "")

ggplot(data = dataVC, aes(x = GrowthForm,y = ksat,
                          fill = as.factor(GrowthForm)))+
  geom_boxplot(alpha = 0.3) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  scale_y_log10() +
  labs(x = "",y= "") +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_blank(),text = element_text(size = 16))

ggsave(plot = last_plot(),
       filename = "./Figures/ksatboxplot.png", dpi = 300, width = 5, height = 5)

# k curves
Top = 35 ; w = 4 ; bottom  = 0.003 ; wbot = 0.0002 ; bottom2 = bottom  - 2*wbot;

bottom3 = 0.1; wbot2 = 0.01
bottom4=0.08; wbot2 = 0.01

yminS05 <- yminS01 <- signif$psi^0*+Top - w
ymaxS05 <- ymaxS01 <- signif$psi^0*+Top + w
yminS05[signif$alpha_k>0.05] <- ymaxS05[signif$alpha_k>0.05] <- NA
ymaxS01[signif$alpha_k>0.01] <- yminS01[signif$alpha_k>0.01] <- NA

ggplot() +
  geom_ribbon(data = bootstrap_sum,aes(x = psi,color = as.factor(GF),
                                       fill = as.factor(GF),ymin = k_low,ymax = k_high),alpha = 0.1,colour = NA) +
  geom_line(data = bootstrap_sum,aes(x = psi,y = k_m,color = as.factor(GF))) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  geom_ribbon(data = data.frame(x = signif$psi,
                                ymin = yminS01,
                                ymax = ymaxS05),aes(x = x,ymin = ymin, ymax = ymax), fill = "darkgrey", colour = NA,alpha = 0.5) +
  geom_ribbon(data = data.frame(x = signif$psi,
                                ymin = yminS05,
                                ymax = ymaxS05),aes(x = x,ymin = ymin, ymax = ymax), fill = "lightgrey", colour = NA,alpha = 0.5) +
  geom_boxploth(data = data.all,aes(x = Pmd,y= 0.03,fill = GrowthForm,color = GrowthForm),alpha = 0.2, outlier.shape = NA,width = 0.4) +
  geom_boxploth(data = data.all,aes(x = Ppd,y= 0.3,fill = GrowthForm,color = GrowthForm),alpha = 0.2, outlier.shape = NA,width = 0.4) +
  geom_boxploth(data = data.all,aes(x = tlp ,y= 0.003,fill = GrowthForm,color = GrowthForm),alpha = 0.2, outlier.shape = NA,width = 0.4) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0.01,0.01)) +
  theme_bw() + theme(legend.position = "none") +
  labs(x = "Water potential",y="Stem conductivity")


data_trait <- data.all %>% dplyr::select(GrowthForm, p50,Pmd,Ppd,tlp) %>% pivot_longer(cols = c(p50,Pmd,Ppd,tlp),names_to = "Trait", values_to = "value") %>% filter(!is.na(value)) %>%
  mutate(Trait = as.factor(Trait))

data_trait$Trait <- factor(data_trait$Trait,levels(data_trait$Trait)[c(3,2,1,4)])

ggplot() +
  geom_ribbon(data = bootstrap_sum,aes(x = psi,color = as.factor(GF),
                                       fill = as.factor(GF),ymin = k_low,ymax = k_high),alpha = 0.1,colour = NA) +
  geom_line(data = bootstrap_sum,aes(x = psi,y = k_m,color = as.factor(GF))) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  geom_ribbon(data = data.frame(x = signif$psi,
                                ymin = yminS01,
                                ymax = ymaxS05),aes(x = x,ymin = ymin, ymax = ymax), fill = "darkgrey", colour = NA,alpha = 0.5) +
  geom_ribbon(data = data.frame(x = signif$psi,
                                ymin = yminS05,
                                ymax = ymaxS05),aes(x = x,ymin = ymin, ymax = ymax), fill = "lightgrey", colour = NA,alpha = 0.5) +
  geom_boxploth(data = data_trait %>% filter(GrowthForm == "Tree"),aes(x = value,y= 0.03, group = Trait),alpha = 0.2, outlier.shape = NA,width = 0.5,
                fill = Cols[2]) +
  geom_boxploth(data = data_trait %>% filter(GrowthForm == "Liana"),aes(x = value,y= 0.005, group = Trait),alpha = 0.2, outlier.shape = NA,width = 0.5,
                fill = Cols[1]) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_log10(expand = c(0.01,0.01)) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 16)) +
  labs(x = "",y="")

ggsave(plot = last_plot(),
       filename = "./Figures/ksatcurve.png", dpi = 300, width = 10, height = 8)


data_trait_fac <- data_trait %>% mutate(t = paste(GrowthForm,Trait,sep = '|'))
data_trait_tree <- data_trait_fac %>% filter(GrowthForm == "Tree")
data_trait_liana <- data_trait_fac %>% filter(GrowthForm == "Liana")

all.pw <- pairwise.wilcox.test(data_trait_fac$value, data_trait_fac$t)$p.value
all.pw[(all.pw < 0.05)] <- NA

tree.pw <- pairwise.wilcox.test(data_trait_tree$value, data_trait_tree$t)$p.value
tree.pw[(tree.pw < 0.01)] <- NA

liana.pw <- pairwise.wilcox.test(data_trait_liana$value, data_trait_liana$t)$p.value
liana.pw[(liana.pw < 0.01)] <- NA

data_trait_fac %>% group_by(t) %>% summarise(v = mean(value,na.rm=TRUE)) %>% arrange(desc(v))

