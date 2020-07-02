rm(list = ls())

library(effsize)
library(dplyr)
library(ggplot2)
library(purrr)
library(scatterpie)
library(LianaHydro)
library(stringr)
library(minpack.lm)

filePV <- "/home/femeunier/Documents/projects/LianaHydro/data/PV.all.csv"
dataPV <- read.csv(filePV,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% filter(Organ == "Leaf") %>% dplyr::select(Species,GrowthForm,p0,tlp,rwc.tlp,epsil,Cft,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Organ) %>%
  mutate(kl=NA,ksat=NA,Al.As=NA,ax=NA,p50=NA,Pmd=NA,Ppd=NA,Id=NA,Function=NA)

fileVC <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
dataVC <- read.csv(fileVC,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% dplyr::select(Species,GrowthForm,kl,ksat,Al.As,ax,p50,Pmd,Ppd,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Id,Function) %>% mutate(
  p0=NA,tlp=NA,rwc.tlp=NA,epsil=NA,Cft=NA,Organ=NA)

data <- rbind(dataPV,dataVC) %>% filter(abs(Lat)<23.5)
data <- correct.obs(data)
data <- add.properties.data(data)

data.filt <- data %>% dplyr::select(GrowthForm,Pmd,Ppd,p50,tlp) %>%
  pivot_longer(c(Pmd,Ppd,p50,tlp), names_to = "Trait", values_to = "value") %>% drop_na(value)

Cols <- c(rgb(0,0,139/255),rgb(0.10,0.50,0.00))

data.fact <-
  data.filt %>% ungroup () %>% mutate(
    GrowthForm_num = case_when(GrowthForm == "Liana" ~ 1,
                               GrowthForm == "Tree"  ~ 2),
    Trait_num = case_when(
      Trait == "Pmd" ~ 2,
      Trait == "Ppd" ~ 4,
      Trait == "p50" ~ 3,
      Trait == "tlp" ~ 1
    )) %>% mutate(pos = case_when(
  GrowthForm_num == 1 ~ ((GrowthForm_num)-1)*4+(Trait_num),
  GrowthForm_num == 2 ~ ((GrowthForm_num))*4-(Trait_num-1)))

ggplot(data = data.fact) +
  geom_boxplot(aes(x = as.factor(pos),y = value,fill = GrowthForm),alpha = 0.5,outlier.colour = NULL) +
  scale_fill_manual(values = Cols) +
  scale_color_manual(values = Cols) +
  labs(y = "Water potential", x = "",fill = "") +
  geom_vline(xintercept = 4.5,linetype = 3) +
  theme_minimal()

data.all <- data.fact %>% mutate(t = as.factor(Trait),
                     gf = as.factor(GrowthForm))

# %>% group_by(t) %>% summarise(p = kruskal.test(formula = value ~ gf)$p.value)

data.liana <- data.all %>% filter(GrowthForm=="Liana")
data.tree <- data.all %>% filter(GrowthForm=="Tree")
pairwise.wilcox.test(data.liana$value, data.liana$t)
pairwise.wilcox.test(data.tree$value, data.tree$t)

all.pw <- pairwise.wilcox.test(data.all$value, data.all$pos)$p.value
all.pw[(all.pw < 0.05)] <- NA

#
# boxplot(formula= Ppd-p12 ~ GrowthForm,data=data)
# summary(aov(formula= Ppd-p12 ~ GrowthForm,data=data))
#
boxplot(formula= Pmd-p50 ~ GrowthForm,data=data)
summary(aov(formula= Pmd-p50 ~ GrowthForm,data=data))

ggplot() +
  geom_point(data = data,aes(x = MAP, y = Pmd-p50,colour = GrowthForm),shape=3)
  geom_point(data = data,aes(x = MAP, y = Ppd-p12,colour = GrowthForm),shape=2)

