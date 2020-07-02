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
library(tidyr)
library(gg3D)

fileVC <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
dataVC <- read.csv(fileVC,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% dplyr::select(Species,GrowthForm,kl,ksat,Al.As,ax,p50,Pmd,Ppd,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Id,Function)

data.all <- correct.obs(dataVC)
data.all <- add.properties.data(data.all)

data.all2 <- data.all %>% mutate(SM50 = Pmd - p50,
                                 SM88 = Pmd - p88) %>% pivot_longer(cols = c(SM50,SM88),
                                        names_to = "SM",
                                        values_to = "SM.value")

ggplot() +
  geom_boxplot(data = data.all2,
               aes(x = Biome,y = SM.value,fill = GrowthForm)) +
  facet_wrap(~ SM,ncol=1) +
  theme_bw()

choat <- read.csv("/home/femeunier/Documents/projects/LianaHydro/data/choatetal2012.csv",stringsAsFactors = FALSE,na.strings=c("","NA"),header=TRUE) %>%
  mutate(Biome = case_when(
    Biome %in% c("tropical rainforest","tropical rain forest") ~ "Tropical rainforest",
    Biome %in% c("tropical seasonal forest") ~ "Tropical seasonal forest",
    TRUE ~ "Other")) %>% filter(Biome != "Other",Growth.form %in% c("Liana","Tree")) %>%
  rename(p50=ψ50,
         p88=ψ88,
         Pmd=ψmin.midday,
         GrowthForm=Growth.form,
         Lat = Latitude,
         Long = Longitude,
         Reference =Reference) %>% dplyr::select(Species,p50,p88,Pmd,GrowthForm,Lat,Long,Biome,Reference) %>% mutate(kl = NA,ksat=NA,Al.As=NA,ax=NA,Ppd=NA,wd=NA,sla=NA,MAP=NA,MAT=NA,Id=NA,Function=NA,p12=NA,
                                                                                                                     p50 = as.numeric(p50),
                                                                                                                     p88 = as.numeric(p88),
                                                                                                                     Pmd = as.numeric(Pmd))

data.all.choat <- rbind(choat,data.all) %>% mutate(SM50 = Pmd - p50,
                                                   SM88 = Pmd - p88)


ggplot() +
  geom_boxplot(data = data.all,
               aes(x = Biome,y = Pmd-p50,fill = GrowthForm)) +
  theme_bw()

ggplot() +
  geom_point(data = data.all,
               aes(x = p50,y = Pmd,color = GrowthForm)) +
  geom_abline(slope=1,intercept = 0,linetype=2) +
  scale_x_continuous(limits = c(-8,0),expand = c(0,0)) +
  scale_y_continuous(limits = c(-8,0),expand = c(0,0)) +
  theme_bw()

# Choat
data.all.choat %>% group_by(GrowthForm,Biome) %>% add_count() %>% summarise(N=mean(n),
                                                                            MAP_m = mean(MAP,na.rm=TRUE))

ggplot() +
  geom_boxplot(data = data.all.choat,
               aes(x = Biome,y = MAP,fill = GrowthForm)) +
  theme_bw()

ggplot() +
  geom_point(data = data.all.choat,
             aes(x = p50,y = Pmd,color = GrowthForm)) +
  geom_abline(slope=1,intercept = 0,linetype=2) +
  # scale_x_continuous(limits = c(-8,0),expand = c(0,0)) +
  # scale_y_continuous(limits = c(-8,0),expand = c(0,0)) +
  theme_bw()


data.SM <- data.all.choat %>% mutate(SM50 = Pmd - p50,
                                     SM88 = Pmd - p88,
                                     SM12 = Ppd - p12) %>% pivot_longer(cols = c(SM12,SM50,SM88),
                                                                        names_to = "SM",
                                                                        values_to = "SM.value")
Cols <- c(rgb(0,0,139/255),rgb(0.10,0.50,0.00))

ggplot(data = data.SM,
       aes(x = GrowthForm,y = SM.value,fill = GrowthForm)) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(position=position_jitterdodge(),alpha=0.3) +
  facet_wrap(~ SM,ncol=3) +
  scale_fill_manual(values = Cols) +
  geom_hline(yintercept = 0,linetype=2) +
  scale_color_manual(values = Cols) +
  labs(x = "",y = "Safety margins [MPa]") +
  theme_light()


data.SM %>% mutate(n = 1) %>% drop_na(SM.value) %>% group_by(GrowthForm,SM) %>% summarise(n = sum(n),
                                                                  p.val = wilcox.test(SM.value)$p.value)
