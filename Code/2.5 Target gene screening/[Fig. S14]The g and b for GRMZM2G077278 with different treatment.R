#==============================================================
#The g and b for GRMZM2G077278 with different treatment
#===============================================================
#Fig.S14

getwd()
source("./Code/Tool/00.load-packages.R")
library(ggpubr)
library(ggsci)

#----
datapath2018<-"./Data/ear_trait_data_2018.rds"
datapath2019<-"./Data/ear_trait_data_2019.rds"


earPhe2018<-read_rds(datapath2018)
earPhe2019<-read_rds(datapath2019)

Phe2018<-earPhe2018%>%
  select(ID1,Genotype)%>%
  unique(.)%>%
  mutate(Genotype=factor(Genotype))
Phe2019<-earPhe2019%>%
  select(ID1,Genotype)%>%
  unique(.)%>%
  mutate(Genotype=factor(Genotype))
rm(earPhe2018)
rm(earPhe2019)

#----
MNormalSeedNum2018 <- readRDS(paste0("./Data/2018fwr-results/","MNormalSeedNum", "_IQR-fwr.rds"))
MNormalSeedNum2019 <- readRDS(paste0("./Data/2019fwr-results/","MNormalSeedNum", "_IQR-fwr.rds"))

g.table2018 <- tibble(Genotype = MNormalSeedNum2018$VARlevels,G = MNormalSeedNum2018$g[, 1])
b.table2018 <- tibble(Genotype = MNormalSeedNum2018$VARlevels,B = MNormalSeedNum2018$b[, 1] + 1)
gb.table2018<-left_join(left_join(g.table2018,b.table2018),Phe2018)

g.table2019 <- tibble(Genotype = MNormalSeedNum2019$VARlevels,G = MNormalSeedNum2019$g[, 1])
b.table2019 <- tibble(Genotype = MNormalSeedNum2019$VARlevels,B = MNormalSeedNum2019$b[, 1] + 1)
gb.table2019<-left_join(left_join(g.table2019,b.table2019),Phe2019)

rm(MNormalSeedNum2018)
rm(MNormalSeedNum2019)
rm(g.table2018)
rm(b.table2018)
rm(g.table2019)
rm(b.table2019)



Test<-gb.table2018%>%
  filter(ID1%in%c("90084","63242"))%>%#"02868"
  mutate(Type=if_else(ID1=="63242","OE","CAS9"))

WT<-gb.table2018%>%
  filter(ID1%in%c("99999"))

p1<-ggplot(data = Test,aes(x=Type,y=G))+
  geom_point(size=2)+
  geom_hline(yintercept = WT$G,linetype=2)+
  labs(x="",y="g")+
  theme_bw()


p2<-ggplot(data = Test,aes(x=Type,y=B))+
  geom_point(size=2)+
  geom_hline(yintercept = WT$B,linetype=2)+
  labs(x="",y="b")+
  theme_bw()

ggpubr::ggarrange(p1,p2)
