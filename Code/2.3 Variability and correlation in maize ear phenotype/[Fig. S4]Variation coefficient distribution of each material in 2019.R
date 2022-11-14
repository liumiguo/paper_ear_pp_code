#=======================================================
#Variation coefficient distribution of each material in 2019
#=======================================================
#Fig. S4

#set path
# setwd(choose.dir())
# getwd()

#Load data
library(tidyverse)

#Read data
list.files("./Data/")
IT2019<-read_rds("./Data/ear_trait_data_2019.rds")
names(IT2019)

#Tidy data
IT2019_cv<-IT2019%>%
  # mutate(ConstructID=str_pad(ConstructID, 5, side = "left", "0"),
  #        EventID=str_pad(EventID, 10, side = "left", "0"))%>%
  within(.,{
    Index2<-NA
    Index2[Index=="MBaldTip"]<-"Bald tip length (mm)"
    Index2[Index=="MBaldTipAreaRatio"]<-"Proportion of bald tip area"
    Index2[Index=="MDiseaseAreaRatio"]<-"Proportion of disease area"
    Index2[Index=="MEarCircumference"]<-"Ear circumference (mm)"
    Index2[Index=="MEarDiameter"]<-"Ear diameter (mm)"
    Index2[Index=="MEmptyAreaRatio"]<-"Proportion of empty area"
    Index2[Index=="MNormalSeedAreaRatio"]<-"Proportion of normal kernel area"
    Index2[Index=="MNormalSeedNum"]<-"Normal kernel number"
    Index2[Index=="MRowSeedNum"]<-"Number of kernels per row"
    Index2[Index=="MSeedThickness"]<-"Kernel thickness (mm)"
    Index2[Index=="MSeedWith"]<-"Kernel width (mm)"
    Index2[Index=="MShrinkageAreaRatio"]<-"Proportion of shrink area"
    Index2[Index=="MShrunkenSeedNum"]<-"Shrunken kernel number"
    Index2[Index=="MSpikeL"]<-"Ear length (mm)"
    Index2[Index=="MSpikeRow"]<-"Number of rows per ear"
    Index2[Index=="MSpikeWidth"]<-"Ear width (mm)"
    Index2[Index=="MSpikeShape"]<-"Ear shape"
  })%>%
  group_by(Index2,Genotype)%>%
  summarise(CV=sd(value)/mean(value))

IT2019_cv_sum<-IT2019_cv%>%
  group_by(Index2)%>%
  summarise(Minv=min(CV),
            Maxv=max(CV),
            JG=Maxv-Minv)

IT2019_cv_WT<-IT2019_cv%>%
  filter(Genotype=="9999901001")
#Plot
ggplot(data = IT2019_cv,aes(x=CV))+
  geom_density(fill="red",color="white",alpha=0.6)+
  geom_vline(data = IT2019_cv_WT,aes(xintercept = CV),color="blue",linetype=2)+
  # geom_histogram()+
  facet_wrap(Index2~.,scales="free")+#free_y
  # scale_y_sqrt()+
  labs(title="2019",x="Coefficient of Variation",y="Density")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())
