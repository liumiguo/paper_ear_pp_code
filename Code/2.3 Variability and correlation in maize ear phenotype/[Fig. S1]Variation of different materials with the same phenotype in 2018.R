#=================================================================
#Variation of different materials with the same phenotype in 2018
#==================================================================
#Distribution map of phenotypic traits of maize ears in 2018
#Fig. S1

#set path
# setwd(choose.dir())
# getwd()

#Load needed package
library(tidyverse)

#Read data
list.files("./Data/")
IT2018<-read_rds("./Data/ear_trait_data_2018.rds")
names(IT2018)

#Tidy data
IT2018_l<-IT2018%>%
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
  filter(!(Index2%in%c("Proportion of disease area","Proportion of shrink area")))


# Batch removal of outliers
# box_out<-function(data){
#   Rowloc<-which(data$value%in%boxplot.stats(data$value)$out)
#   todel<-sort(unique(unlist(Rowloc)))
#   re<-data[-todel,]
#   return(re)
# }
# 
# IT2018_l2<-plyr::ddply(.data = IT2018_l,.variables = c("Index2"),.fun = box_out)


#plot
unique(IT2018_l$Index2)
names(IT2018_l)

ggplot(data = IT2018_l,aes(x=value))+
  # geom_density(fill="red",color="white",alpha=0.6)+
  geom_histogram()+
  facet_wrap(Index2~.,scales="free")+
  # scale_y_sqrt()+
  labs(x=" ",y="Count")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())

#Save data after removing outliers
# write_rds(IT2018_l2,"./Data/2018_ear_trait_long_data_rmNA.rds")

