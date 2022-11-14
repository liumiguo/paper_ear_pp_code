#=============================================
#Correlation analysis of ear phenotype in 2019
#==============================================
#Fig. S5

#set path
# setwd(choose.dir())
# getwd()

#Load needed package
library(tidyverse)

#Read data
list.files("./Data/")
ear2019<-read_rds("./Data/ear_trait_data_2019.rds")
names(ear2019)

#Tidy data
ear2019_IT_index<-ear2019%>%
  spread(.,key = Index,value = value)%>%
  drop_na()


names(ear2019_IT_index)[5:22]<-c("Rep","ID2","IndexRep","Bald tip length",
                                 "Proportion of bald tip area",
                                 "Ear circumference",
                                 "Ear diameter",
                                 "Proportion of empty area",
                                 "Proportion of normal kernel area",
                                 "Normal kernel number",
                                 "Number of kernels per row",
                                 "Kernel thickness",
                                 "Kernel width",
                                 "Shrunken kernel number",
                                 "Ear length",
                                 "Number of rows per ear",
                                 "Ear shape",
                                 "Ear width")
ear2019_IT_index$ID<-1:nrow(ear2019_IT_index)
ear2019_IT_index_l<-ear2019_IT_index%>%
  gather(.,key = Index,value = value,`Bald tip length`:`Ear width`)

# #Batch removal of outliers
# box_out<-function(data){
#   Rowloc<-which(data$value%in%boxplot.stats(data$value)$out)
#   todel<-sort(unique(unlist(Rowloc)))
#   value2<-data$value
#   value2[todel]<-NA
#   value2[is.na(value2)]=mean(value2,na.rm=T)
#   data$value2<-value2
#   # re<-data[-todel,]
#   re<-data
#   return(re)
# }
# 
# names(ear2019_IT_index_l)
# ear2019_IT_index_l2<-plyr::ddply(.data = ear2019_IT_index_l,.variables = c("Index"),.fun = box_out)

#Tidy data after removing outliers
names(ear2019_IT_index_l)
ear2019_IT_index_spr<-ear2019_IT_index_l%>%
  spread(.,key = Index,value = value)


#plot
library(ggcorrplot)
names(ear2019_IT_index_spr)
corr2019 <- round(cor(ear2019_IT_index_spr[,c(9:23)]), 2)
head(corr2019[, 1:6])

p2019 <- cor_pmat(ear2019_IT_index_spr[,c(9:23)])
head(p2019[, 1:4])

ggcorrplot(corr2019)
ggcorrplot(corr2019,outline.color = "white")
ggcorrplot(corr2019,type = "lower", outline.color = "white")#下三角形
ggcorrplot(corr2019,type = "lower", lab = TRUE)
ggcorrplot(corr2019,type = "full", p.mat = p2019)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text.x = element_text(size=14,angle = 45, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"))

