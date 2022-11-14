#========================================================================
#2018-2019 number of excellent g and excellent b genes
#========================================================================
#Fig 3A

# setwd("./Data/")
getwd()
# getwd()

library(tidyverse)
library(VennDiagram)
library(data.table)
source("./Code/Tool/plot.FW2.R")

#file path
datapath<-"./Data/ear_trait_data_2018.rds"

list.files("./Data/2018fwr-results")

#input data
NormalSeed2018 <- readRDS("./Data/2018fwr-results/MNormalSeedNum_IQR-fwr.rds")
Earlength2018 <- readRDS("./Data/2018fwr-results/MSpikeL_IQR-fwr.rds")
RowSeedNum2018 <- readRDS("./Data/2018fwr-results/MRowSeedNum_IQR-fwr.rds")
# NormalSeed2019 <- readRDS("./2019fwr-results/MNormalSeedNum_IQR-fwr.rds")

#
Phe2018<-read_rds(datapath)%>%
  select(ID1,Genotype)%>%
  unique(.)%>%
  mutate(Genotype=factor(Genotype))
# Phe2019<-read_csv(datapath2019)%>%
#   select(CAUID,ConstructID,EventID)%>%
#   unique(.)%>%
#   rename(Genotype=EventID)


##normal kernel number
g.MNormalSeedNum2018 <- tibble(Genotype = NormalSeed2018$VARlevels,G = NormalSeed2018$g[, 1])
b.MNormalSeedNum2018 <- tibble(Genotype = NormalSeed2018$VARlevels,B = NormalSeed2018$b[, 1] + 1)
# a_b<-left_join(g.MNormalSeedNum2018,b.MNormalSeedNum2018)
gb.MNormalSeedNum2018<-left_join(left_join(g.MNormalSeedNum2018,b.MNormalSeedNum2018),Phe2018)

##ear length
g.Earlength2018 <- tibble(Genotype = Earlength2018$VARlevels,G = Earlength2018$g[, 1])
b.Earlength2018 <- tibble(Genotype = Earlength2018$VARlevels,B = Earlength2018$b[, 1] + 1)
gb.Earlength2018<-left_join(left_join(g.Earlength2018,b.Earlength2018),Phe2018)

##number of kernel per row
g.RowSeedNum2018 <- tibble(Genotype = RowSeedNum2018$VARlevels,G = RowSeedNum2018$g[, 1])
b.RowSeedNum2018 <- tibble(Genotype = RowSeedNum2018$VARlevels,B = RowSeedNum2018$b[, 1] + 1)
gb.RowSeedNum2018<-left_join(left_join(g.RowSeedNum2018,b.RowSeedNum2018),Phe2018)

# rm(Phe2018)
# rm(Phe2019)

#The g value of screened strains was higher than that of WT strains

gb.MNormalSeedNum2018_num<-gb.MNormalSeedNum2018%>%
  select(Genotype,ID1,G,B)%>%
  mutate(WT_g=4.43910156,#
         WT_b=1.0979779)%>%
  filter(ID1!="99999")%>%
  mutate(betterG=G>WT_g)%>%
  filter(betterG==TRUE)%>%# 
  mutate(lowerB=B<=WT_b,
         upB=B>WT_b)%>%
  group_by(ID1)%>%
  summarise(NSN_GbetterNum=sum(betterG),
            # NSN_GTotalNum=n(),
            NSN_BlowerNum=sum(lowerB),
            NSN_upNum=sum(upB))%>%
  arrange(-NSN_GbetterNum)#

names(gb.MNormalSeedNum2018_num)
gb.MNormalSeedNum2018_num_sum<-gb.MNormalSeedNum2018_num%>%
  group_by(NSN_GbetterNum,NSN_upNum)%>%#
  summarise(N=n())

p2018<-ggplot(data = gb.MNormalSeedNum2018_num_sum,aes(x=NSN_GbetterNum,y=NSN_upNum))+
  geom_point(aes(size=N),color="blue")+
  geom_text(aes(x=NSN_GbetterNum+0.2,label=N),size=5)+
  geom_vline(xintercept = 2.5,linetype=2,color="red")+
  labs(x="Number of events with\n g value greater than WT\n",
       y="Number of events with\n b value greater than WT\n",
       size=" The number of\n genes:")+
  theme_bw()+
  theme(axis.text.x = element_text(size=14,angle = 0, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"))


#-----2019
datapath2019<-"./Data/ear_trait_data_2019.rds"

list.files("./Data/2019fwr-results")

#
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


g.Earlength2019 <- tibble(Genotype = Earlength2019$VARlevels,G = Earlength2019$g[, 1])
b.Earlength2019 <- tibble(Genotype = Earlength2019$VARlevels,B = Earlength2019$b[, 1] + 1)
gb.Earlength2019<-left_join(left_join(g.Earlength2019,b.Earlength2019),Phe2019)


g.RowSeedNum2019 <- tibble(Genotype = RowSeedNum2019$VARlevels,G = RowSeedNum2019$g[, 1])
b.RowSeedNum2019 <- tibble(Genotype = RowSeedNum2019$VARlevels,B = RowSeedNum2019$b[, 1] + 1)
gb.RowSeedNum2019<-left_join(left_join(g.RowSeedNum2019,b.RowSeedNum2019),Phe2019)

# rm(Phe2018)
# rm(Phe2019)


gb.MNormalSeedNum2019_num<-gb.MNormalSeedNum2019%>%
  select(Genotype,ID1,G,B)%>%
  mutate(WT_g=3.955502042,
         WT_b=1.18289164)%>%
  filter(ID1!="99999")%>%
  mutate(betterG=G>WT_g)%>%
  filter(betterG==TRUE)%>%
  mutate(lowerB=B<=WT_b,
         upB=B>WT_b)%>%
  group_by(ID1)%>%
  summarise(NSN_GbetterNum=sum(betterG),
            # NSN_GTotalNum=n(),
            NSN_BlowerNum=sum(lowerB),
            NSN_upNum=sum(upB))%>%
  arrange(-NSN_GbetterNum)

names(gb.MNormalSeedNum2019_num)
gb.MNormalSeedNum2019_num_sum<-gb.MNormalSeedNum2019_num%>%
  group_by(NSN_GbetterNum,NSN_upNum)%>%
  summarise(N=n())
p2019<-ggplot(data = gb.MNormalSeedNum2019_num_sum,aes(x=NSN_GbetterNum,y=NSN_upNum))+
  geom_point(aes(size=N),color="blue")+
  geom_text(aes(x=NSN_GbetterNum+0.25,label=N),size=5)+
  geom_vline(xintercept = 2.5,linetype=2,color="red")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6))+
  scale_y_continuous(breaks = c(1,2,3,4,5,6))+
  labs(x="Number of events with\n g value greater than WT\n",
       y="Number of events with\n b value greater than WT\n",
       size=" The number of\n genes:")+
  theme_bw()+
  theme(axis.text.x = element_text(size=14,angle = 0, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"))


library(ggpubr)
ggarrange(p2018,p2019,nrow = 1,common.legend = TRUE,labels = c("2018","2019"),
          hjust = -4,
          vjust = 3,
          font.label = list(size = 10))
