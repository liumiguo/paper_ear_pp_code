#=================================================
#2018
#=================================================
library(dplyr)
library(agricolae)
library(readr)
library(tidyr)
library(stringr)
# library(data.table)
# Load and prep data 
list.files()
earPhe<-read_rds("./Data/ear_trait_data_2018.rds")
earPhe2019<-earPhe%>%
  mutate(Environment=str_c(Site,"_",Rep))%>%
  filter(Index=="MNormalSeedNum")%>%
  select(Genotype,Environment,IndexRep,value)%>%
  drop_na()

rm(earPhe)
# earPhe2019_2<-earPhe2019%>%
#   group_by(EventID,Environment,IndexRep)%>%
#   summarise(M=mean(value))


#----AMMI model
names(earPhe2019)
str(earPhe2019)
gc()

memory.limit(size = 9999999999999) 
model<- AMMI(ENV=earPhe2019$Environment,
             GEN=earPhe2019$EventID,
             REP=1,
             Y=earPhe2019$value,console=F,PC=F)
write_rds(model,"2018AMMI_reslut.rds")
##查看方差分析
model$ANOVA

##查看GEI部分的主成分
model$analysis#查看主成分显著不显著

#查看ASV和YSI结果：
#ASV AMMI stability value
#YSI Yield stability index
#means average genotype by environment
Idx2<-index.AMMI(model)

write.csv(Idx2,"2018AMMI_index.csv")
