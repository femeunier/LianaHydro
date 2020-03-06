rm(list = ls())

library(ggplot2)

file <- "/home/femeunier/Documents/projects/LianaHydro/data/VC.all.csv"
data <- read.csv(file)


mLianas <- nlsLM(data = data %>% filter(GrowthForm == "Liana"),
            ax ~ a*(-p50)**b,
            start=list(a=54.4, b=-1.17), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                               printEval = TRUE, warnOnly = TRUE))
mTrees <- nlsLM(data = data %>% filter(GrowthForm == "Tree"),
                 ax ~ a*(-p50)**b,
                 start=list(a=54.4, b=-1.17), control = nls.control(maxiter = 500, tol = 1e-05, minFactor = 1/1024/10,
                                                                    printEval = TRUE, warnOnly = TRUE))

data_added <- data %>% mutate(ax.mod = case_when(
  GrowthForm == "Liana" & is.na(ax) ~  coef(mLianas)[1]*(-p50)** coef(mLianas)[2],
  GrowthForm == "Tree" & is.na(ax) ~  coef(mTrees)[1]*(-p50)** coef(mTrees)[2],
  TRUE ~ ax))

psi = seq(-10,0,length.out = 100)
colors <- data.frame(growthform = c("Liana","Tree"),
                     col = c("Darkblue","darkgreen"))

plot(NA,NA,xlim=c(min(psi),max(psi)),ylim=c(0,100))

data_uni <- data_added %>% filter(!is.na(p50),!is.na(ax.mod))
lianas <- trees <- c()

for (i in seq(1,nrow(data_uni))){
  p50 <- data_uni$p50[i]
  ax50 <- data_uni$ax.mod[i]
  ksat <- data_uni$ksat[i]
  GF <-  data_uni$GrowthForm[i]
  if (!is.na(p50) & !is.na(ax50)){

    b = p50
    a = -4*ax50*p50/100
    PLC = polynomial2(psi,a,b)
    lines(psi,PLC,col = as.character(colors$col[match(GF,colors$growthform)]))

    # b = -p50
    # a = 4*ax50/100
    # PLC = weibull(psi,b,a)
    # lines(psi,PLC,col = as.character(colors$col[match(GF,colors$growthform)]))

    if (GF == "Liana"){
      lianas <- rbind(lianas,PLC)
    } else {
      trees <- rbind(trees,PLC)
    }
  }
}

Nbootstrap = 250
Nliana <- nrow(lianas)
Ntree <- nrow(trees)

ksat_liana <- data_added %>% filter(!is.na(ksat) & GrowthForm == "Liana")%>%pull(ksat)
Nliana_K <- length(ksat_liana)

ksat_tree <- data_added %>% filter(!is.na(ksat) & GrowthForm == "Tree")%>%pull(ksat)
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
                     function.names = c("weibull","sigmoidal","polynomial","polynomial2"))
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

signif <- bootstrap %>% group_by(psi) %>% summarise(PLC = kruskal.test(formula = PLC ~ GF)$p.value,
                                                    k  = kruskal.test(formula = k ~ GF)$p.value)

bootstrap_sum <- bootstrap %>% filter(k > 0)  %>% group_by(GF,psi) %>% summarise(PLC_m = mean(PLC),
                                                              PLC_low = quantile(PLC,0.025),
                                                              PLC_high = quantile(PLC,0.975),
                                                              k_m = mean(k),
                                                              k_low = quantile(k,0.025),
                                                              k_high = quantile(k,0.975))

Cols <- c(rgb(0,0,139/255),rgb(0.10,0.50,0.00))

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
  # scale_y_log10(expand = c(0.01,0.01)) +
  theme_bw()

# P50 boxplot
ggplot(data = data.summary, aes(x = GF,y = P50, fill = as.factor(GF)))+
  geom_boxplot(alpha = 0.2) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  theme_bw()

ggplot(data = data_added, aes(x = GrowthForm,y = p50,
                              fill = as.factor(GrowthForm)))+
  geom_boxplot(alpha = 0.3) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  labs(x = "") +
  theme_bw() + theme(legend.position = "none")


# ksat boxplot
ggplot(data = data_added, aes(x = GrowthForm,y = ksat,
                              fill = as.factor(GrowthForm)))+
  geom_boxplot(alpha = 0.3) +
  scale_color_manual(values = Cols) +
  scale_fill_manual(values = Cols) +
  scale_y_log10() +
  labs(x = "") +
  theme_bw() + theme(legend.position = "none")

# Pmd boxplot
quantile(data_added %>% filter(GrowthForm=="Tree" & !is.na(Pmd)) %>% pull(Pmd),c(0.025,0.975))
quantile(data_added %>% filter(GrowthForm=="Liana" & !is.na(Pmd)) %>% pull(Pmd),c(0.025,0.975))

########################################################################################
#
# plot(psi,colMeans(trees),type='l',col='darkgreen')
# polygon(c(psi,rev(psi)),
#         c(colMeans(trees)-1.96*apply(trees,2,sd)/sqrt(nrow(trees)),
#           rev(colMeans(trees)+1.96*apply(trees,2,sd)/sqrt(nrow(trees)))),
#         col=rgb(0,1,0,0.1),border = NA)
#
# lines(psi,colMeans(lianas),type='l',col='darkblue')
# polygon(c(psi,rev(psi)),
#         c(colMeans(lianas)-1.96*apply(lianas,2,sd)/sqrt(nrow(lianas)),
#           rev(colMeans(lianas)+1.96*apply(lianas,2,sd)/sqrt(nrow(lianas)))),
#         col=rgb(0,0,1,0.1),border = NA)
#
########################################################################################
#
# plot(psi,colMeans(treesK),type='l',col='darkgreen',log='y',ylim=c(0.1,20))
# polygon(c(psi,rev(psi)),
#         c(colMeans(treesK)-1.96*apply(treesK,2,sd)/sqrt(nrow(treesK)),
#           rev(colMeans(treesK)+1.96*apply(treesK,2,sd)/sqrt(nrow(treesK)))),
#         col=rgb(0,1,0,0.1),border = NA)
#
#
# lines(psi,colMeans(lianasK),type='l',col='darkblue')
# polygon(c(psi,rev(psi)),
#         c(colMeans(lianasK)-1.96*apply(lianasK,2,sd)/sqrt(nrow(lianasK)),
#           rev(colMeans(lianasK)+1.96*apply(lianasK,2,sd)/sqrt(nrow(lianasK)))),
#         col=rgb(0,0,1,0.1),border = NA)
