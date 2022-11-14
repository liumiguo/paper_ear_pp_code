#===========================================
#the accuracy for under full ear
#===========================================

getwd()

library(tidyverse)

list.files("./Code/2.7 Accuracy/")

data<-read_csv("./Data/plumped ears.csv")

names(data)
accuracy<-data%>%
  filter(MV!=0)%>%
  mutate(acc=100-(abs((MV-RV))/RV)*100)%>%
  group_by(Index)%>%
  summarise(Accm=mean(acc))


