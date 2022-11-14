#===============================================================================
#200个检测果穗分组的精度
#===============================================================================

getwd()
# setwd("./Code/2.7 Accuracy/")

library(tidyverse)
library(data.table)

list.files()

data<-fread("./Code/2.7 Accuracy/ear_traits_group.csv")


ear_pre_measure_group_acc<-data%>%
  drop_na()%>%
  filter(MV!=0)%>%
  mutate(acc=100-(abs((MV-RV))/RV)*100)%>%
  group_by(Index,Group)%>%
  summarise(Accm=mean(acc))%>%
  spread(.,key = Group,value = Accm)


