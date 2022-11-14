#=====================================================
#Global prediction performance
#=====================================================
#Table 1

getwd()

library(tidyverse)
library(readxl)

xlsxfile="./Data/Data_Acc.xlsx"
readxl::excel_sheets(xlsxfile)


#Accuracy
data_acc<-read_xlsx(xlsxfile,"Accuracy")%>%
  drop_na()%>%
  mutate(Acc=100-abs(RV-MV)/RV*100)%>%
  group_by(Index)%>%
  summarise(Accuracy=mean(Acc))

#Trueness
data_tru<-read_xlsx(xlsxfile,"Trueness")%>%
  filter(RV!=0)%>%
  drop_na()%>%
  mutate(MVm=(MV1+MV2+MV3)/3)%>%
  mutate(Acc=100-abs(RV-MVm)/RV*100)%>%
  group_by(Index)%>%
  summarise(Trueness=mean(Acc))

#Precision
data_pre<-read_xlsx(xlsxfile,"Precision")%>%
  filter(RV!=0)%>%
  drop_na()%>%
  gather(.,key = Time,value = value,MV1:MV3)%>%
  group_by(Index,ID)%>%
  summarise(SD=sd(value),
            M=mean((value)))%>%
  mutate(Pre=100-SD/M*100)%>%
  group_by(Index)%>%
  summarise(Precision=mean(Pre))

#Reproducibility
data_Repro<-read_xlsx(xlsxfile,"Reproducibility")%>%
  drop_na()%>%
  gather(.,key = Time,value = value,MV1:MV5)%>%
  group_by(Index,ID)%>%
  summarise(M=mean(value),
            SD=sd(value))%>%
  mutate(CV=SD/M)%>%
  group_by(Index)%>%
  summarise(Reproducibility=100-mean(CV)*100)

#Stability
data_Stab<-read_xlsx(xlsxfile,"Stability")%>%
  drop_na()%>%
  gather(.,key = Days,value = value,Day1:Day10)%>%
  group_by(Index,ID)%>%
  summarise(M=mean(value),
            SD=sd(value))%>%
  mutate(CV=SD/M)%>%
  group_by(Index)%>%
  summarise(Stability=100-mean(CV)*100)

#"Technical repeatability"
data_Tec<-read_xlsx(xlsxfile,"Technical repeatability")%>%
  drop_na()%>%
  gather(.,key = Index,value = value,3:8)%>%
  filter(value!=0)%>%
  mutate(loc=str_locate(ID,"_")[,1],
         ID1=str_sub(ID,1,loc-1))%>%
  drop_na()%>%
  group_by(Index,ID1)%>%
  summarise(M=mean(value),
            SD=sd(value))%>%
  mutate(CV=SD/M)%>%
  group_by(Index)%>%
  summarise(TechnicalRepeatability=100-mean(CV)*100)



#Join data
join_data<-Reduce(left_join,list(data_acc,
                  data_tru,
                  data_pre,
                  data_Repro,
                  data_Stab,
                  data_Tec))

# write_csv(join_data,"./Code/2.7 Accuracy/Table1.csv")

