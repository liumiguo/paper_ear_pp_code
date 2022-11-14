#==========================================================
#EXCELLENT GENES IN 2019
#===========================================================
# setwd("./Data/")
getwd()

library(tidyverse)
library(VennDiagram)
source("./Code/Tool/plot.FW2.R")

datapath2019<-"Data/ear_trait_data_2019.rds"

list.files("./Data/2019fwr-results")


NormalSeed2019 <- readRDS("./Data/2019fwr-results/MNormalSeedNum_IQR-fwr.rds")
Earlength2019 <- readRDS("./Data/2019fwr-results/MSpikeL_IQR-fwr.rds")
RowSeedNum2019 <- readRDS("./Data/2019fwr-results/MRowSeedNum_IQR-fwr.rds")
# NormalSeed2019 <- readRDS("./2019fwr-results/MNormalSeedNum_IQR-fwr.rds")


Phe2019<-read_rds(datapath2019)%>%
  select(ID1,Genotype)%>%
  unique(.)%>%
  mutate(Genotype=factor(Genotype))
# Phe2019<-read_csv(datapath2019)%>%
#   select(CAUID,ConstructID,EventID)%>%
#   unique(.)%>%
#   rename(Genotype=EventID)


g.MNormalSeedNum2019 <- tibble(Genotype = NormalSeed2019$VARlevels,G = NormalSeed2019$g[, 1])
b.MNormalSeedNum2019 <- tibble(Genotype = NormalSeed2019$VARlevels,B = NormalSeed2019$b[, 1] + 1)
gb.MNormalSeedNum2019<-left_join(left_join(g.MNormalSeedNum2019,b.MNormalSeedNum2019),Phe2019)


gb.MNormalSeedNum2019_betterG<-gb.MNormalSeedNum2019%>%
  select(Genotype,ID1,G,B)%>%
  mutate(WT_g=3.955502042,
         WT_b=1.18289164)%>%
  # filter(ConstructID!="99999")%>%
  mutate(betterG=G>WT_g)%>%
  filter(betterG==TRUE)


gb.MNormalSeedNum2019_betterG_num<-gb.MNormalSeedNum2019_betterG%>%
  mutate(lowerB=B<=WT_b,
         upB=B>WT_b)%>%
  group_by(ID1)%>%
  summarise(NSN_GbetterNum=sum(betterG),
            # NSN_GTotalNum=n(),
            NSN_BlowerNum=sum(lowerB),
            NSN_upNum=sum(upB))%>%
  arrange(-NSN_GbetterNum)%>%
  filter(NSN_GbetterNum>=3)




Construct_ID<-"93009"

par(mfrow=c(1,3))
NormalSeed2019_Genotype<-gb.MNormalSeedNum2019_betterG%>%
  filter(ID1%in%c(Construct_ID,"99999"))
plot.FW2(x=NormalSeed2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal seed \n number in 2019"))
plot.FW2(x=Earlength2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2019"))
plot.FW2(x=RowSeedNum2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number kernels\n per row in 2019"))

par(mfrow=c(1,3),mar=c(4,4,4,4))
plot.new()
NormalSeed2019_Genotype<-gb.MNormalSeedNum2019%>%
  filter(ID1%in%c(Construct_ID,"99999"))
plot.FW2(x=NormalSeed2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal seed \n number in 2019"))
plot.FW2(x=Earlength2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2019"))
plot.FW2(x=RowSeedNum2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number kernels\n per row in 2019"))






# ConstructID<-unique(gb.MNormalSeedNum2019_betterG_num$ConstructID)
ConstructID<-c("00725",
               "02241",
               "02445",
               "02479",
               "02926",
               "03522",
               "03773",
               "03796",
               "03803",
               "03813",
               "80017",
               "80024",
               "80027",
               "80035",
               "80051",
               "80067",
               "80073",
               "80080")

for (i in 1:length(ConstructID)) {
  # i=1
  Construct_ID<-ConstructID[i]
  png(str_c("./Code/2.5 Target gene screening/2019gene/",Construct_ID,".png"),width=816,height=336)
  par(mfrow=c(1,3))
  NormalSeed2019_Genotype<-gb.MNormalSeedNum2019_betterG%>%
    filter(ID1%in%c(Construct_ID,"99999"))
  plot.FW2(x=NormalSeed2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal kernel \n number in 2019"))
  plot.FW2(x=Earlength2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2019"))
  plot.FW2(x=RowSeedNum2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number of kernels\n per row in 2019"))
  dev.off()
}
