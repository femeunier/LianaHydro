rm(list = ls())

library(effsize)
library(dplyr)
library(ggplot2)
library(purrr)
library(scatterpie)

file <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
data <- read.csv(file,stringsAsFactors = FALSE,na.strings=c("","NA")) %>% mutate(kl = log10(kl),
                                                                                 ksat = log10(ksat))

vars <- c("kl","ksat","Al.As","ax","p50","Pmd","Ppd","wd","sla","MAP","MAT")

Cohensd <- data.frame()

rmin = 0.05
rmax = 0.15

results <-
  map_dfr(vars,function(var){
    data.temp <- data %>% select(Species,GrowthForm,var,Reference) %>% distinct()

    GF <- data.temp %>% filter(!is.na(!!as.symbol(var)) &
                               is.finite(!!as.symbol(var))) %>% group_by(GrowthForm) %>% add_count() %>% summarise(N = mean(n))

    d <- effsize::cohen.d(formula = as.formula(paste(var," ~ GrowthForm")),data = data.temp,na.rm=TRUE,noncentral=TRUE)

    KWtest <- kruskal.test(data.temp %>%pull(!!as.symbol(var)),
                              as.factor(data.temp %>%pull(GrowthForm)))$p.value

    return(data.frame(variable = as.factor(var),
                        cohensd = as.numeric(d$estimate),
                        cohenlow = as.numeric(d$conf.int[1]),
             cohenhigh = as.numeric(d$conf.int[2]),
             pvalue = KWtest,
             nliana = GF %>% filter(GrowthForm == "Liana") %>% pull(N),
             ntree = GF %>% filter(GrowthForm == "Tree") %>% pull(N)))
      }) %>% mutate(N = nliana + ntree) %>% add_rownames("id") %>% mutate(id = as.numeric(id),
                                                                          r = rmin + rmax*(N-min(N))/(max(N) - min(N)),
                                                                          signif = case_when(pvalue < 0.001 ~ "***",
                                                                                             pvalue < 0.01 ~ "**",
                                                                                             pvalue < 0.05 ~ "*",
                                                                                             TRUE ~ "N.S."))

Nmin = min(results$N)
Nmax = max(results$N)

ggplot(data = results) +
  geom_errorbarh(data = results,aes(xmin = cohenlow,xmax = cohenhigh,
                                    y = -id),height = 0.) +
  geom_scatterpie(aes(y=-id, x=cohensd,r = r),data = results,
                  cols=c("nliana", "ntree")) +
  scale_fill_manual(values = c("darkblue","darkgreen")) +
  geom_vline(xintercept = 0,linetype=2) +
  scale_y_continuous(breaks = seq(-1,-length(vars),-1),labels = vars,name = "",sec.axis = dup_axis(labels = results$signif)) +
  scale_x_continuous(name = "Cohen d [-]") +
  geom_scatterpie_legend((results$r), x = -1, y = -1,n = 2,
                         labeller = function(x) x= round(1/100*(Nmin + (Nmax-Nmin)*(x- rmin)/(rmax-rmin)))*100) +
  theme_bw() + theme(legend.position = "none")


Tref <- unique(data %>% filter(!is.na(p50) & GrowthForm=="Tree") %>% pull(Reference))
Lref <- unique(data %>% filter(!is.na(p50) & GrowthForm=="Liana") %>% pull(Reference))
I <- intersect(Tref,Lref)

effsize::cohen.d(formula = p50 ~ GrowthForm,data= data %>% filter(Reference == I[5]),
                 na.rm=TRUE,paired = FALSE,hedges.correction = TRUE,noncentral=TRUE)
