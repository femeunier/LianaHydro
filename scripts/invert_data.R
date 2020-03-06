rm(list=ls())

library(dplyr)
library(LianaHydro)
library(minpack.lm)

data.file <- file.path(getwd(),"data","Acacia.csv")
data <- read.csv(data.file,stringsAsFactors=FALSE,
                 na.strings=c("","NA")) %>% dplyr::select(Datum,Time,WP..MPa.,height..m.,time.0.1ml..s.,length..m.,Surface.in..m..,Surface.out..m..) %>%
  rename(date = Datum,
         t = Time,
         wp = WP..MPa.,
         h = height..m.,
         tQ = time.0.1ml..s.,
         l = length..m.,
         Sin = Surface.in..m..,
         Sout = Surface.out..m..)

dataf <- fill.NA.data(data) %>% mutate(din = sqrt(Sin/1000/pi)*2*1000*2/2,
                                       dout = sqrt(Sout/1000/pi)*2*1000*2/2,
                                       P = h*9.81*1000,
                                       tQ = as.numeric(tQ),
                                       Q = 0.1/1000/1000/tQ)

datat <- dataf %>% group_by(date,t) %>% summarise(K = coef(lm(Q ~ P))[2],
                                                  r2 = summary(lm(Q ~ P))$adj.r.squared,
                                                  din = mean(din),
                                                  Sinprim = mean(pi*(din**2)/4)/1/1e6,
                                                  dout = mean(dout),
                                                  Soutprim = mean(pi*(dout**2)/4)/1e6,
                                                  lm = mean(l),
                                                  k = 1000*K/((Sinprim+Soutprim)/2)*lm,
                                                  wp_m = mean(wp)) %>% ungroup() %>%
  mutate(PLC = 100*(1 - k/k[1])) %>% arrange(desc(wp_m))  %>% rename(psi = wp_m)


# plot((datat$din+datat$dout)/2,datat$k)
# plot(datat$psi,datat$PLC)
# temp <- datat %>% select(Sinprim,K,k,PLC,psi)
#
# k_DBH <- function(psi,Sinprim,param1,param2,a,b){
#   Ks = a + b*Sinprim
#   PLC = weibull(psi,param1,param2)
#   K = (1-PLC/100)*Ks
#   return(K)
# }
#
# m <- nlsLM(data = temp,
#            K ~ k_DBH(psi,Sinprim,param1,param2,a,b),
#            start=list(param1 =  2,
#                       param2 = 2,
#                       a = 0,
#                       b = 0.001),
#            lower = c(0,0.01,0,0),
#            upper = c(Inf,Inf,Inf,Inf),
#            control = nls.control(maxiter = 1000, tol = 1e-10, minFactor = 1/1024/10,
#                                  printEval = TRUE, warnOnly = TRUE))
# #
# #
# Coef <- coef(m)
# plot(temp$Sinprim,k_DBH(temp$psi,temp$Sinprim,Coef[1],Coef[2],Coef[3],Coef[4]))
# plot(temp$psi,polynomial(temp$psi,Coef[1],Coef[2]))
#
# plot(temp$psi,temp$K,pch=19)
# lines(temp$psi,k_DBH(temp$psi,temp$din,Coef[1],Coef[2],Coef[3],Coef[4]),pch=19,type='p',col="red")
# lines(temp$psi,(1-best.model[[1]]$PLC.predict/100)*(temp$K[1]),pch=19,type='p',col="blue")
#
# plot(temp$psi,temp$PLC,pch=19,ylim=c(0,100))
# lines(temp$psi,polynomial(temp$psi,Coef[1],Coef[2]),pch=19,type='p',col="red")
# lines(temp$psi,(best.model[[1]]$PLC.predict),pch=19,type='p',col="blue")

data2opt <- datat %>% dplyr::select(psi,PLC)

All.models <-
  opt.data(data2opt,
           function.names = c("weibull","sigmoidal","polynomial","polynomial2","cumnorm"))
All.models <- add.properties(All.models,x = c(12,50,88))
best.model <-
  find.best.model(All.models)

plot(data2opt$psi,data2opt$PLC,ylim=c(0,max(max(data2opt$PLC),100)),
     pch=19,xlab = "Psi",ylab = "PLC")
plot.models(All.models,col=c('black'),add=TRUE,highlight=TRUE)

