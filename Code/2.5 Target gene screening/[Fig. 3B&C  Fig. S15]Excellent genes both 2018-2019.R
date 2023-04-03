#=====================================================================================
#10 genes shared in 2018-2019 (including Wayne map)
#=====================================================================================
#Fig. 3B&C

# setwd("./Data/")
getwd()
# getwd()

library(tidyverse)
library(VennDiagram)
library(data.table)
source("./Code/Tool/plot.FW3.R")


datapath<-"./Data/ear_trait_data_2018.rds"
datapath2019<-"./Data/ear_trait_data_2019.rds"

list.files("./Data/2018fwr-results")

#-----2018

NormalSeed2018 <- readRDS("./Data/2018fwr-results/MNormalSeedNum_IQR-fwr.rds")
Earlength2018 <- readRDS("./Data/2018fwr-results/MSpikeL_IQR-fwr.rds")
RowSeedNum2018 <- readRDS("./Data/2018fwr-results/MRowSeedNum_IQR-fwr.rds")

Phe2018<-read_rds(datapath)%>%
  select(ID1,Genotype)%>%
  unique(.)%>%
  mutate(Genotype=factor(Genotype))

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

#----2019
NormalSeed2019 <- readRDS("./Data/2019fwr-results/MNormalSeedNum_IQR-fwr.rds")
Earlength2019 <- readRDS("./Data/2019fwr-results/MSpikeL_IQR-fwr.rds")
RowSeedNum2019 <- readRDS("./Data/2019fwr-results/MRowSeedNum_IQR-fwr.rds")
# NormalSeed2019 <- readRDS("./2019fwr-results/MNormalSeedNum_IQR-fwr.rds")


Phe2019<-read_rds(datapath2019)%>%
  select(ID1,Genotype)%>%
  unique(.)%>%
  mutate(Genotype=factor(Genotype))

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
nrow(gb.MNormalSeedNum2019_betterG_num)

G_list<-list("Candidate genes\n from 2018"=gb.MNormalSeedNum2018_betterG_num$ID1,
             "Candidate genes\n from 2019"=gb.MNormalSeedNum2019_betterG_num$ID1)

library(ggvenn)
ggvenn(G_list,c("Candidate genes\n from 2018","Candidate genes\n from 2019"),
       # fill_color = c("blue", "red"),
       fill_color = c("blue", "red"),
       fill_alpha = 0.6,
       set_name_size = 8,
       stroke_color = "grey",
       stroke_alpha=0.5,
       stroke_size = 0.5,
       text_size=7,
       text_color = "white")+
  # labs(title = "The number of excellent genes")+
  theme(plot.title = element_text(hjust = 0.5,size=20))



#---The FWR plot of 10 genes 
Construct_ID<-"00344"

par(mfrow=c(2,3))
NormalSeed2018_Genotype<-gb.MNormalSeedNum2018_betterG%>%
  filter(ID1%in%c(Construct_ID,"99999"))
plot.FW2(x=NormalSeed2018,plotVAR = c(NormalSeed2018_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal seed \n number in 2018"))
plot.FW2(x=Earlength2018,plotVAR = c(NormalSeed2018_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2018"))
plot.FW2(x=RowSeedNum2018,plotVAR = c(NormalSeed2018_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number kernels\n per row in 2018"))

NormalSeed2019_Genotype<-gb.MNormalSeedNum2019_betterG%>%
  filter(ID1%in%c(Construct_ID,"99999"))
plot.FW2(x=NormalSeed2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal seed \n number in 2019"))
plot.FW2(x=Earlength2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2019"))
plot.FW2(x=RowSeedNum2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number kernels\n per row in 2019"))

ConstructID<-c("00344","01299","01924","02454","02900","02916","03824","80072","90084","91078")
for (i in 1:length(ConstructID)) {
  Construct_ID<-ConstructID[i]
  png(str_c("./Code/2.5 Target gene screening/10gene/",Construct_ID,".png"),width=816,height=672)
  par(mfrow=c(2,3))
  NormalSeed2018_Genotype<-gb.MNormalSeedNum2018_betterG%>%
    filter(ID1%in%c(Construct_ID,"99999"))
  plot.FW2(x=NormalSeed2018,plotVAR = c(NormalSeed2018_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal kernel \n number in 2018"))
  plot.FW2(x=Earlength2018,plotVAR = c(NormalSeed2018_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2018"))
  plot.FW2(x=RowSeedNum2018,plotVAR = c(NormalSeed2018_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number of kernels\n per row in 2018"))

  NormalSeed2019_Genotype<-gb.MNormalSeedNum2019_betterG%>%
    filter(ID1%in%c(Construct_ID,"99999"))
  plot.FW2(x=NormalSeed2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Normal kernel \n number in 2019"))
  plot.FW2(x=Earlength2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Ear length in 2019"))
  plot.FW2(x=RowSeedNum2019,plotVAR = c(NormalSeed2019_Genotype$Genotype,9999901001),main=str_c(Construct_ID,": Number of kernels\n per row in 2019"))

  dev.off()
}



