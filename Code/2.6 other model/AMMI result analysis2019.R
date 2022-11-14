#=============================================================
#AMMI result analysis in 2019
#=============================================================

# setwd(choose.dir())

library(tidyverse)
library(paletteer)

list.files("./Code/2.6 other model/")

#---2019
Index2019<-read_csv("./Code/2.6 other model/result/2019AMMI_index.csv")
EarTrait2019<-read_rds("./Data/ear_trait_data_2019.rds")
EarTrait2019<-EarTrait2019%>%
  filter(Index=="MNormalSeedNum")%>%
  group_by(ID1,Genotype)%>%
  summarise(CV=sd(value)/mean(value))

names(Index2019)[1]<-"Genotype"

Join2019<-left_join(Index2019,EarTrait2019)


Select2019<-c("00344","01299","01924","02454","02900","02916","03824","80072","90084","91078",
              "00725","02241","02445","02479","02926","03522","03773","03796","03803","03813",
              "80017","80024","80027","80035","80051","80067","80073","80080")

Join2019_select<-Join2019%>%
  filter(ID1%in%Select2019)

#means
WT_means<-Join2019%>%filter(Genotype=="9999901001")%>%.$means
Join2019_select_means<-Join2019_select%>%
  filter(means>WT_means)%>%
  mutate(Rate=(means-WT_means)/WT_means*100)


str(Join2019_select_means)
p_mean<-ggplot(data = Join2019_select_means,aes(x=ID1,y=means,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_means,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="Mean value of normal kernel number",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30))

#CV
WT_CV<-Join2019%>%filter(Genotype=="9999901001")%>%.$CV
Join2019_select_CV<-Join2019_select_means%>%
  mutate(Rate=(CV-WT_CV)/WT_CV*100)
p_cv<-ggplot(data = Join2019_select_CV,aes(x=ID1,y=CV,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_CV,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="CV of normal kernel number",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30))

#ASV
WT_ASV<-Join2019%>%filter(Genotype=="9999901001")%>%.$ASV
Join2019_select_ASV<-Join2019_select_means%>%
  mutate(Rate=(ASV-WT_ASV)/WT_ASV*100)
p_ASV<-ggplot(data = Join2019_select_ASV,aes(x=ID1,y=ASV,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_ASV,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="ASV\n",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x= element_text(angle = 30))

#YSI
WT_YSI<-Join2019%>%filter(Genotype=="9999901001")%>%.$YSI
Join2019_select_YSI<-Join2019_select_means%>%
  mutate(Rate=(YSI-WT_YSI)/WT_YSI*100)
p_YSI<-ggplot(data = Join2019_select_YSI,aes(x=ID1,y=YSI,color=ID1))+
  geom_point(size=3)+
  geom_hline(yintercept = WT_YSI,linetype=2)+
  # scale_color_manual(values=paletteer_d("dichromat::SteppedSequential_5"))+
  labs(x="",y="YSI",title="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30))

ggpubr::ggarrange(p_mean,p_cv,p_ASV,ncol=1)
