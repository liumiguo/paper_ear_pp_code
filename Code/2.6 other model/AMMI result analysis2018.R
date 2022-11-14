#=============================================================
#AMMI result analysis in 2018
#=============================================================

# setwd(choose.dir())

library(tidyverse)
library(paletteer)

list.files("./Code/2.6 other model/")

#---2018
Index2018<-read_csv("./Code/2.6 other model/result/2018AMMI_index.csv")
EarTrait2018<-read_rds("./Data/ear_trait_data_2018.rds")
names(EarTrait2018)
EarTrait2018<-EarTrait2018%>%
  filter(Index=="MNormalSeedNum")%>%
  group_by(ID1,Genotype)%>%
  summarise(CV=sd(value)/mean(value))

names(Index2018)[1]<-"Genotype"

Join2018<-left_join(Index2018,EarTrait2018)

Select2018<-c("00344","01299","01924","02454","02900","02916","03824","80072","90084","91078",
              "01627","01835","01875","02157","02840","91047")

Join2018_select<-Join2018%>%
  filter(ID1%in%Select2018)

#means
WT_means<-Join2018%>%filter(Genotype=="9999901001")%>%.$means
Join2018_select_means<-Join2018_select%>%
  filter(means>WT_means)%>%
  mutate(Rate=(means-WT_means)/WT_means*100)


str(Join2018_select_means)
p_mean<-ggplot(data = Join2018_select_means,aes(x=ID1,y=means,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_means,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="Mean value of normal kernel number",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30))

#CV
WT_CV<-Join2018%>%filter(Genotype=="9999901001")%>%.$CV
Join2018_select_CV<-Join2018_select_means%>%
  mutate(Rate=(CV-WT_CV)/WT_CV*100)
p_cv<-ggplot(data = Join2018_select_CV,aes(x=ID1,y=CV,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_CV,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="CV of normal kernel number",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30))

#ASV
WT_ASV<-Join2018%>%filter(Genotype=="9999901001")%>%.$ASV
Join2018_select_ASV<-Join2018_select_means%>%
  mutate(Rate=(ASV-WT_ASV)/WT_ASV*100)
p_ASV<-ggplot(data = Join2018_select_ASV,aes(x=ID1,y=ASV,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_ASV,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="ASV\n",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x= element_text(angle = 30))

#YSI
WT_YSI<-Join2018%>%filter(Genotype=="9999901001")%>%.$YSI
Join2018_select_YSI<-Join2018_select_means%>%
  mutate(Rate=(YSI-WT_YSI)/WT_YSI*100)
p_YSI<-ggplot(data = Join2018_select_YSI,aes(x=ID1,y=YSI,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_YSI,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="YSI",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30))

ggpubr::ggarrange(p_mean,p_cv,p_ASV,ncol=1)
