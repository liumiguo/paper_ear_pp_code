#===========================================================
#2018-2019 QCD merge plot
#===========================================================
# setwd("./Data/")
getwd()
library(tidyverse)

list.files(".")


p2018<-readRDS("./Data/2018_QCD2.rds")
p2019<-readRDS("./Data/2019_QCD2.rds")


ggpubr::ggarrange(p2018,p2019,common.legend = T,labels = c("2018","2019"),
                  hjust = -2,
                  vjust = 2.5)
