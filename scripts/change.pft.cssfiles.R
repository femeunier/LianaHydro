rm(list = ls())

library(dplyr)
library(ggplot2)

directory <- "/home/femeunier/Downloads/nph14009-sup-0003-notess2/PV_demography"
directory.out <- "/home/femeunier/Downloads/nph14009-sup-0003-notess2/PV_demography_pft"
if(!dir.exists(directory.out)) {dir.create(directory.out)}

site <- "PV_2009_2P_site6.lat10.4lon-85.4"
site.out <- "PV_2009_2P_site6_liana.lat10.4lon-85.4"

css.file <- file.path(directory,paste0(site,".css"))
data.cohorts <- read.csv(css.file,header = TRUE,stringsAsFactors = FALSE,sep = " ")

data.cohorts <- data.cohorts %>% mutate(pft = case_when(
  pft == 24 ~ 2,
  pft == 26 ~ 2,
  pft == 27 ~ 3,
  pft == 29 ~ 4))

# we could add some lianas here
data.cohorts_liana <- rbind(c(0,1,9998,5,35,17,0.1,0,0,0),
                            data.cohorts,
                            c(0,2,9999,3,33,17,0.1,0,0,0))


ggplot(data.cohorts_liana) +
  geom_histogram(aes(x = dbh,fill = as.factor(pft),colour = as.factor(pft)),position = 'stack')

write.table(data.cohorts,file = file.path(directory.out,paste0(site,".css")),sep = " ",row.names = FALSE)
system2("cp",paste(file.path(directory,paste0(site,".pss")),file.path(directory.out,paste0(site,".pss"))))

write.table(data.cohorts_liana,file = file.path(directory.out,paste0(site.out,".css")),sep = " ",row.names = FALSE)
system2("cp",paste(file.path(directory,paste0(site,".pss")),file.path(directory.out,paste0(site.out,".pss"))))

