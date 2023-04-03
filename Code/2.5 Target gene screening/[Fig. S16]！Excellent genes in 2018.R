#==========================================================
#Excellent genes in 2018
#===========================================================
# setwd("./Data/")
getwd()
# getwd()

library(tidyverse)
library(VennDiagram)
source("./Code/Tool/plot.FW3.R")


datapath<-"Data/ear_trait_data_2018.rds"

list.files("./Data/2018fwr-results")


NormalSeed2018 <- readRDS("./Data/2018fwr-results/MNormalSeedNum_IQR-fwr.rds")
Earlength2018 <- readRDS("./Data/2018fwr-results/MSpikeL_IQR-fwr.rds")
RowSeedNum2018 <- readRDS("./Data/2018fwr-results/MRowSeedNum_IQR-fwr.rds")
# NormalSeed2019 <- readRDS("./2019fwr-results/MNormalSeedNum_IQR-fwr.rds")


Phe2018<-read_rds(datapath)%>%
  select(ID1,Genotype)%>%
  unique(.)%>%
  mutate(Genotype=factor(Genotype))
# Phe2019<-read_csv(datapath2019)%>%
#   select(CAUID,ConstructID,EventID)%>%
#   unique(.)%>%
#   rename(Genotype=EventID)


g.MNormalSeedNum2018 <- tibble(Genotype = NormalSeed2018$VARlevels,G = NormalSeed2018$g[, 1])
b.MNormalSeedNum2018 <- tibble(Genotype = NormalSeed2018$VARlevels,B = NormalSeed2018$b[, 1] + 1)
gb.MNormalSeedNum2018<-left_join(left_join(g.MNormalSeedNum2018,b.MNormalSeedNum2018),Phe2018)


gb.MNormalSeedNum2018_betterG<-gb.MNormalSeedNum2018%>%
  select(Genotype,ID1,G,B)%>%
  mutate(WT_g=4.43910156,
         WT_b=1.0979779)%>%
  # filter(ConstructID!="99999")%>%
  mutate(betterG=G>=WT_g)%>%
  filter(betterG==TRUE)


gb.MNormalSeedNum2018_betterG_num<-gb.MNormalSeedNum2018_betterG%>%
  mutate(lowerB=B<=WT_b,
         upB=B>WT_b)%>%
  group_by(ID1)%>%
  summarise(NSN_GbetterNum=sum(betterG),
            # NSN_GTotalNum=n(),
            NSN_BlowerNum=sum(lowerB),
            NSN_upNum=sum(upB))%>%
  arrange(-NSN_GbetterNum)%>%
  filter(NSN_GbetterNum>=3)

# write_csv(gb.MNormalSeedNum2018_betterG_num,"./2018年最终筛选的优良基因.csv")


Construct_ID<-"00344"

par(mfrow=c(1,3))
NormalSeed2019_Genotype<-gb.MNormalSeedNum2018_betterG%>%
  filter(ID1%in%c(Construct_ID,"99999"))
plot.FW2(x=NormalSeed2018,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal seed \n number in 2018"))
plot.FW2(x=Earlength2018,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2018"))
plot.FW2(x=RowSeedNum2018,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number kernels\n per row in 2018"))

ConstructID<-unique(gb.MNormalSeedNum2018_betterG_num$ID1)

ConstructID<-c("01627",
               "01835",
               "01875",
               "02157",
               "02840",
               "91047")
for (i in 1:length(ConstructID)) {
  Construct_ID<-ConstructID[i]
  png(str_c("./Code/2.5 Target gene screening/2018gene/",Construct_ID,".png"),width=816,height=336)
  par(mfrow=c(1,3))
  NormalSeed2019_Genotype<-gb.MNormalSeedNum2018_betterG%>%
    filter(ID1%in%c(Construct_ID,"99999"))
  plot.FW2(x=NormalSeed2018,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal kernel \n number in 2018"))
  plot.FW2(x=Earlength2018,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2018"))
  plot.FW2(x=RowSeedNum2018,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number of kernels\n per row in 2018"))
  dev.off()
}
