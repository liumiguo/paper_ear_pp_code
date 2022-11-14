#===========================================================================
#34 excellent genes
#===========================================================================
#Fig. 4


getwd()
source("./Code/Tool/00.load-packages.R")
library(ggpubr)
library(ggsci)
library(data.table)
library(readr)

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
join_2019<-left_join(g.table2019,b.table2019)
gb.table2019<-left_join(join_2019,Phe2019)

rm(MNormalSeedNum2018)
rm(MNormalSeedNum2019)
rm(g.table2018)
rm(b.table2018)
rm(g.table2019)
rm(b.table2019)
#----34 genes
G34<-read_csv("./Data/34genes_info.csv")

names(G34)
G34_form<-G34%>%
  mutate(MaterialID=sprintf('%05d',MaterialID))

G10<-G34_form%>%
  filter(Year=="2018-2019")

#----10 genes

gb.table2018_wt<-gb.table2018%>%
  filter(Genotype=="9999901001")
gb.table2019_wt<-gb.table2019%>%
  filter(Genotype=="9999901001")


gb.table2018_beter<-gb.table2018%>%
  filter(ID1%in%unique(G10$MaterialID)[1:10])%>%
  mutate(G_frac=(G-gb.table2018_wt$G)/gb.table2018_wt$G*100,
         B_frac=(B-gb.table2018_wt$B)/gb.table2018_wt$B*100)%>%
  filter(G_frac>=0)%>%
  arrange(ID1,B_frac,G_frac)

gb.table2019_beter<-gb.table2019%>%
  filter(ID1%in%unique(G10$MaterialID)[1:10])%>%
  mutate(G_frac=(G-gb.table2019_wt$G)/gb.table2019_wt$G*100,
         B_frac=(B-gb.table2019_wt$B)/gb.table2019_wt$B*100)%>%
  filter(G_frac>=0)%>%
  arrange(ID1,B_frac,G_frac)


#---
gb.table2018_beter_l<-gb.table2018_beter%>%
  select(ID1,G_frac,B_frac)%>%
  gather(.,key = Index,value = frac,G_frac:B_frac)

gb.table2019_beter_l<-gb.table2019_beter%>%
  select(ID1,G_frac,B_frac)%>%
  gather(.,key = Index,value = frac,G_frac:B_frac)

gb.table2018_beter_l$Year="2018"
gb.table2019_beter_l$Year="2019"

gb_better_2018_2019_l<-rbind(gb.table2018_beter_l,gb.table2019_beter_l)%>%
  mutate(Index=if_else(Index=="G_frac","g","b"))

p2y<-ggplot(data = gb_better_2018_2019_l, aes(x=frac,y=Index))+
  geom_point(aes(shape=Index,color=Year),size=4,alpha=0.6)+
  geom_vline(xintercept=0,linetype=3)+
  # scale_color_npg()+
  scale_color_manual(values=c("red","blue"))+
  facet_grid(ID1~.)+
  labs(x="\nIncreased ratio relative to wild type (%)",y="")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none")+
  theme(axis.text.x = element_text(size=14,angle = 0, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"
        ), 
        strip.text.y = element_text(
          size = 12, color = "black", face = "bold.italic"
        ) )


#-----
G24<-G34_form%>%
  filter(Year!="2018-2019")

#----
##WT
gb.table2018_wt<-gb.table2018%>%
  filter(Genotype=="9999901001")
gb.table2019_wt<-gb.table2019%>%
  filter(Genotype=="9999901001")

##INCREACTION FRACTION
gb.table2018_beter<-gb.table2018%>%
  filter(ID1%in%c(G24%>%filter(Year=="2018")%>%.$MaterialID))%>%
  mutate(G_frac=(G-gb.table2018_wt$G)/gb.table2018_wt$G*100,
         B_frac=(B-gb.table2018_wt$B)/gb.table2018_wt$B*100)%>%
  filter(G_frac>=0)%>%
  arrange(ID1,B_frac,G_frac)

gb.table2019_beter<-gb.table2019%>%
  filter(ID1%in%c(G24%>%filter(Year=="2019")%>%.$MaterialID))%>%
  mutate(G_frac=(G-gb.table2019_wt$G)/gb.table2019_wt$G*100,
         B_frac=(B-gb.table2019_wt$B)/gb.table2019_wt$B*100)%>%
  filter(G_frac>=0)%>%
  arrange(ID1,B_frac,G_frac)



#---
gb.table2018_beter_l<-gb.table2018_beter%>%
  select(ID1,G_frac,B_frac)%>%
  gather(.,key = Index,value = frac,G_frac:B_frac)

gb.table2019_beter_l<-gb.table2019_beter%>%
  select(ID1,G_frac,B_frac)%>%
  gather(.,key = Index,value = frac,G_frac:B_frac)

gb.table2018_beter_l$Year="2018"
gb.table2019_beter_l$Year="2019"

gb_better_2018_2019_l<-rbind(gb.table2018_beter_l,gb.table2019_beter_l)%>%
  mutate(Index=if_else(Index=="G_frac","g","b"))

##2018
p3<-ggplot(data = gb_better_2018_2019_l%>%filter(Year =="2018"),aes(x=frac,y=Index))+
  geom_point(aes(shape=Index,color=Year),size=3,alpha=0.6)+
  geom_vline(xintercept=0,linetype=2)+
  # scale_color_npg()+
  scale_color_manual(values=c("red","blue"))+
  # scale_x_continuous(breaks = c(-50,0,100,200,300,400))+
  xlim(-50,400)+
  facet_grid(ID1~.)+
  labs(x="\n",y="")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none")+
  theme(axis.text.x = element_text(size=14,angle = 0, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"
        ), 
        strip.text.y = element_text(
          size = 12, color = "black", face = "bold.italic"
        ) )

##2019
l2019<-gb_better_2018_2019_l%>%filter(Year =="2019")


p1<-ggplot(data = gb_better_2018_2019_l%>%filter(Year =="2019",ID1%in%unique(l2019$ID1)[1:9]),aes(x=frac,y=Index))+
  geom_point(aes(shape=Index,color=Year),size=3,alpha=0.6)+
  geom_vline(xintercept=0,linetype=2)+
  # scale_color_npg()+
  scale_color_manual(values=c("blue"))+
  xlim(-50,400)+
  facet_grid(ID1~.,)+
  # labs(x="\nIncreased ratio relative to wild type (%)",y="")+
  labs(x="\n",y="")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none")+
  theme(axis.text.x = element_text(size=14,angle = 0, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"
        ), 
        strip.text.y = element_text(
          size = 12, color = "black", face = "bold.italic"
        ) )

p2<-ggplot(data = gb_better_2018_2019_l%>%filter(Year =="2019",ID1%in%unique(l2019$ID1)[10:18]),aes(x=frac,y=Index))+
  geom_point(aes(shape=Index,color=Year),size=3,alpha=0.6)+
  geom_vline(xintercept=0,linetype=2)+
  # scale_color_npg()+
  scale_color_manual(values=c("blue"))+
  xlim(-50,400)+
  facet_grid(ID1~.,)+
  labs(x="\n",y="")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none")+
  theme(axis.text.x = element_text(size=14,angle = 0, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"),
        strip.text.x = element_text(
          size = 12, color = "black", face = "bold.italic"
        ),
        strip.text.y = element_text(
          size = 12, color = "black", face = "bold.italic"
        ) )

ggpubr::ggarrange(p2y,p3,p1,p2,ncol = 4)






