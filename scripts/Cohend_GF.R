rm(list = ls())

library(effsize)
library(dplyr)
library(ggplot2)
library(purrr)
library(scatterpie)
library(LianaHydro)
library(minpack.lm)
library(stringr)

filePV <- "/home/femeunier/Documents/projects/LianaHydro/data/PV.all.csv"
dataPV <- read.csv(filePV,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% filter(Organ == "Leaf") %>% dplyr::select(Species,GrowthForm,p0,tlp,rwc.tlp,epsil,Cft,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Organ) %>%
  mutate(kl=NA,ksat=NA,Al.As=NA,ax=NA,p50=NA,Pmd=NA,Ppd=NA,Id=NA,Function=NA)

fileVC <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
dataVC <- read.csv(fileVC,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% dplyr::select(Species,GrowthForm,kl,ksat,Al.As,ax,p50,Pmd,Ppd,wd,sla,MAP,MAT,Reference,Biome,Long,Lat,Id,Function) %>% mutate(
  p0=NA,tlp=NA,rwc.tlp=NA,epsil=NA,Cft=NA,Organ=NA)

data <- rbind(dataPV,dataVC) %>% filter(abs(Lat)<23.5)
data <- correct.obs(data)
data <- add.properties.data(data)

vars <- c("kl","ksat","Al.As","ax","p50","p12","p88","Pmd","Ppd","tlp","wd","sla","MAP","MAT")
hlines <- -(c(3,7,10,12)+0.5)

Cohensd <- data.frame()

rmin = 0.05
rmax = 0.15

results <-
  map_dfr(vars,function(var){
    data.temp <- data %>% select(Species,GrowthForm,var,Reference) %>% distinct() %>% filter(!is.na(!!as.symbol(var)) &
                                                                                               is.finite(!!as.symbol(var)))
    Nstudyliana <- length(unique(data.temp %>% filter(GrowthForm == "Liana") %>% pull(Reference)))
    Nstudytree <- length(unique(data.temp %>% filter(GrowthForm == "Tree") %>% pull(Reference)))

    GF <- data.temp %>% group_by(GrowthForm) %>% add_count() %>% summarise(N = mean(n))

    d <- effsize::cohen.d(formula = as.formula(paste(var," ~ GrowthForm")),data = data.temp,na.rm=TRUE,noncentral=TRUE)

    KWtest <- kruskal.test(data.temp %>%pull(!!as.symbol(var)),
                           as.factor(data.temp %>%pull(GrowthForm)))$p.value

    return(data.frame(variable = as.factor(var),
                      cohensd = as.numeric(d$estimate),
                      cohenlow = as.numeric(d$conf.int[1]),
                      cohenhigh = as.numeric(d$conf.int[2]),
                      pvalue = as.numeric(KWtest),
                      nliana = as.numeric(GF %>% filter(GrowthForm == "Liana") %>% pull(N)),
                      ntree = as.numeric(GF %>% filter(GrowthForm == "Tree") %>% pull(N)),
                      Nstudyliana = as.numeric(Nstudyliana),
                      Nstudytree = as.numeric(Nstudytree)))
  }) %>% mutate(N = nliana + ntree,
                Nstudy = Nstudyliana + Nstudytree) %>% add_rownames("id") %>% mutate(id = as.numeric(id),
                                                                                     r = rmin + rmax*(N-min(N))/(max(N) - min(N)),
                                                                                     signif = case_when(pvalue < 0.001 ~ "***",
                                                                                                        pvalue < 0.01 ~ "**",
                                                                                                        pvalue < 0.05 ~ "*",
                                                                                                        TRUE ~ "N.S."))

Nmin = min(results$N)
Nmax = max(results$N)

Cols <- c(rgb(0,0,139/255),rgb(0.10,0.50,0.00))
ggplot(data = results) +
  geom_errorbarh(data = results,aes(xmin = cohenlow,xmax = cohenhigh,
                                    y = -id),height = 0.) +
  geom_scatterpie(aes(y=-id, x=cohensd,r = r),data = results,
                  cols=c("nliana", "ntree")) +
  scale_fill_manual(values = Cols) +
  geom_vline(xintercept = 0,linetype=3) +
  geom_hline(yintercept = hlines,linetype=2) +
  scale_y_continuous(breaks = seq(-1,-length(vars),-1),labels = vars,name = "",sec.axis = dup_axis(labels = results$signif)) +
  scale_x_continuous(name = "Cohen d [-]") +
  geom_scatterpie_legend((results$r), x = -1, y = -1,n = 2,
                         labeller = function(x) x= round(1/100*(Nmin + (Nmax-Nmin)*(x- rmin)/(rmax-rmin)))*100) +
  theme_bw() + theme(legend.position = "none",
                     text = element_text(size = 16))

ggsave(plot = last_plot(),
       filename = "./Figures/Cohensd.png", dpi = 300, width = 6, height = 10)


# Tref <- unique(data %>% filter(!is.na(p50) & GrowthForm=="Tree") %>% pull(Reference))
# Lref <- unique(data %>% filter(!is.na(p50) & GrowthForm=="Liana") %>% pull(Reference))
# I <- intersect(Tref,Lref)
#
# effsize::cohen.d(formula = p50 ~ GrowthForm,data= data %>% filter(Reference == I[5]),
#                  na.rm=TRUE,paired = FALSE,hedges.correction = TRUE,noncentral=TRUE)
