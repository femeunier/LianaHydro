rm(list = ls())

library(dplyr)
library(ggplot2)
library(purrr)
library(LianaHydro)
library(tidyr)
library(scatterpie)
library(stringr)
library(RColorBrewer)

Latmin <- 90

filePV <- "/home/femeunier/Documents/projects/LianaHydro/data/PV.all.csv"
dataPV <- read.csv(filePV,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% filter(Organ == "Leaf") %>% dplyr::select(Species,GrowthForm,p0,tlp,rwc.tlp,epsil,Cft,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Organ) %>%
  mutate(kl=NA,ksat=NA,Al.As=NA,ax=NA,p50=NA,Pmd=NA,Ppd=NA,Id=NA,Function=NA)

fileVC <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
dataVC <- read.csv(fileVC,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% dplyr::select(Species,GrowthForm,kl,ksat,Al.As,ax,p50,Pmd,Ppd,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Id,Function) %>% mutate(
  p0=NA,tlp=NA,rwc.tlp=NA,epsil=NA,Cft=NA,Organ=NA)
data <- rbind(dataVC,dataPV)

struct <- c("Species","GrowthForm","MAP","MAT","Reference","Biome","Long","Lat")
vars <- c("p0","tlp","kl","ksat","Al.As","ax","p50","Pmd","Ppd")

traits <- data %>% dplyr::select(c(struct,vars))

traits_temp <- traits %>% select(one_of(vars))
traits_temp[is.na(traits_temp)] <- 0
traits_temp[traits_temp!=0] <- 1

traits_fin <- cbind(traits %>% select(one_of(struct)),
                    traits_temp)

traits.species <- traits_fin %>% gather(trait, present,vars)

worldmap <- map_data("world")

rmin = 1
rmax = 6

trait.summ <- traits.species %>% filter(abs(Lat)<Latmin) %>% group_by(Long,Lat,GrowthForm) %>% summarise(tot.trait = sum(present))
trait.summ.GF <- trait.summ %>% spread(GrowthForm,tot.trait) %>% replace_na(list(Liana = 0,Tree = 0)) %>% mutate(N = Liana + Tree) %>% ungroup() %>%
  mutate(r = rmin + rmax*(N-min(N))/(max(N)-min(N)))

data <- trait.summ.GF
data <- merge.sources(data,tol=2)

Cols <- c(rgb(0,0,139/255),rgb(0.10,0.50,0.00))
Nmin = min(data$N)
Nmax = max(data$N)

ggplot() +
  geom_map(data = worldmap, map = worldmap, aes(x=long, y=lat, map_id=region), col = "lightgrey", fill = "lightgrey") +
  geom_scatterpie(data = data,
                  aes(x=Long, y=Lat, r = r),cols=c("Liana","Tree"),alpha = 0.8,color = NA) +
  scale_fill_manual(values = Cols) +
  scale_y_continuous(limits = c(-1,1)*(Latmin+2.5)) +
  scale_x_continuous(limits = c(-1,1)*180,expand=c(0,0)) +
  labs(x="",y="",fill="") +
  geom_scatterpie_legend((data$r), x = -150, y = -20,n = 3,
                         labeller = function(x) x= pmax(1,round(1/50*(Nmin + (Nmax-Nmin)*(x- rmin)/(rmax-rmin)))*50)) +
  theme_minimal() + theme(panel.grid = element_blank(),
                          axis.text = element_blank(),
                          legend.position = c(0.12,0.3),
                          text = element_text(size=14))

ggsave(plot = last_plot(),
      filename = "./Figures/map.png", dpi = 999, width = 10, height = 5)

traits_g.species <- traits.species %>% mutate(trait_g = as.factor(case_when(
  trait %in% c("kl","ksat","Al.As") ~ "Conductivity",
  trait %in% c("ax","p50") ~ "VC",
  trait %in% c("Pmd","Ppd","tlp","p0") ~ "LWP",
  trait %in% c("wd","sla") ~ "structural",
  TRUE ~ "Others")))

traits_g.species$trait_g <- factor(traits_g.species$trait_g,levels(traits_g.species$trait_g)[c(1,4,2,3)])

Observations <- traits_g.species %>% group_by(GrowthForm,trait_g) %>% summarise(N = sum(present))

Observations <- Observations %>% ungroup() %>% mutate(trait_fac = as.factor(trait_g),
                                                      GF_fac = as.factor(GrowthForm))
levels(Observations$trait_fac) = seq(1,length(levels(Observations$trait_fac)))
levels(Observations$GF_fac) = c(1,2)
Observations <- Observations %>% mutate(palette = (as.numeric(GF_fac)-1)*3 + as.numeric(trait_fac))

trtCol <- (brewer.pal(length(levels(Observations$trait_fac)), "Blues"))
conCol <- (brewer.pal(length(levels(Observations$trait_fac)), "Greens"))
Palette<-c(trtCol,conCol)

ggplot(Observations, aes(GrowthForm, N, fill = trait_g))+
  geom_bar(position = "stack",stat = "identity",fill = Palette) +
  coord_flip() + theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(x="",y="Number of trait.species")

ggsave(plot = last_plot(),
       filename = "./Figures/Ntraits.png", dpi = 300, width = 20, height = 5)

# Number of studies and species
Nsum <- traits.species %>% filter(present != 0) %>% group_by(GrowthForm) %>% summarise(study = length(unique(Reference)),
                                                                                       species = length(unique(Species)),
                                                                                       trait = sum(present)) %>%
  pivot_longer(c(study,species,trait), names_to = "type", values_to = "N") %>% mutate(type = as.factor(type),
                                                                                      GrowthForm = as.factor(GrowthForm))

Nsum$type <- factor(Nsum$type,levels(Nsum$type)[c(2,1,3)])
Nsum$GrowthForm <- factor(Nsum$GrowthForm,levels(Nsum$GrowthForm)[c(2,1)])

ggplot(Nsum, aes(GrowthForm, N, fill = GrowthForm))+
  geom_bar(position = "identity",stat = "identity",show.legend = FALSE) +
  facet_wrap(~ type,scales= "free_x",ncol=1) + theme_minimal() +
  scale_fill_manual(values = rev(Cols)) +
  coord_flip() +labs(x = "",fill="",y="") + theme(strip.text = element_blank(),
                                                  text = element_text(16))

ggsave(plot = last_plot(),
       filename = "./Figures/Numbers.png", dpi = 300, width = 5, height = 2)

# MAP
Observations_MAP <- traits_g.species %>% mutate(MAP_g = case_when(
  MAP < 1500 ~ 1,
  MAP > 2500 ~ 3,
  TRUE ~ 2
)) %>% filter(present != 0) %>% group_by(GrowthForm,MAP_g,trait_g) %>% summarise(N = sum(present))

ggplot(Observations_MAP, aes(GrowthForm, N, fill = trait_g))+
  geom_bar(position = "stack",stat = "identity",fill = rep(Palette,3)) +
  coord_flip() + theme_minimal() +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank()) +
  facet_wrap(~MAP_g,nrow = 3) +
  labs(x="",y="") + theme(strip.text = element_blank(),
                          text = element_text(16))

ggsave(plot = last_plot(),
       filename = "./Figures/MAP_N.png", dpi = 300, width = 5, height = 2)
