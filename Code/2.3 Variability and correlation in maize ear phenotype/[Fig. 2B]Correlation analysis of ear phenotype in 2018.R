#=============================================
#Correlation analysis of ear phenotype in 2018
#==============================================
#A heat map was used to show the correlation between ear phenotypic traits
#Fig. 2B

#Set working path
# setwd(choose.dir())
# getwd()

#Load needed package
library(tidyverse)
library(PerformanceAnalytics)

#Read data
list.files("./Data/")
ear2018<-read_rds("./Data/ear_trait_data_2018.rds")
names(ear2018)


#Tidy data
ear2018_IT_index<-ear2018%>%
  spread(.,key = Index,value = value)%>%
  drop_na()


names(ear2018_IT_index)[5:22]<-c("Rep","ID2","IndexRep","Bald tip length",
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
                     "Ear width"
                     )
ear2018_IT_index$ID<-1:nrow(ear2018_IT_index)
ear2018_IT_index_l<-ear2018_IT_index%>%
  gather(.,key = Index,value = value,`Bald tip length`:`Ear width`)

#Remove outlier value
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
# names(ear2018_IT_index_l)
# ear2018_IT_index_l2<-plyr::ddply(.data = ear2018_IT_index_l,.variables = c("Index"),.fun = box_out)

#Tidy Data after removing outliers
ear2018_IT_index_l2<-ear2018_IT_index_l
names(ear2018_IT_index_l2)
ear2018_IT_index_spr<-ear2018_IT_index_l2%>%
  mutate(Group=str_c(Site,MatClass,ID1,Genotype,Rep,sep = "_"))%>%
  spread(.,key = Index,value = value)

names(ear2018_IT_index_spr)

#Calculate correlation matrix
library(ggcorrplot)
names(ear2018_IT_index_spr)
corr2018 <- round(cor(ear2018_IT_index_spr[,c(10:24)]), 2)
head(corr2018[, 1:6])

p2018 <- cor_pmat(ear2018_IT_index_spr[,c(10:24)])
head(p2018[, 1:4])

ggcorrplot(corr2018)
ggcorrplot(corr2018,outline.color = "white")
ggcorrplot(corr2018,type = "lower", outline.color = "white")
ggcorrplot(corr2018,type = "lower", lab = TRUE)


mytheme <- theme(axis.ticks.length.y = unit(-0.15,"cm"), 
                 axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm")),
                 axis.line = element_line(size = 0.8),
                 axis.ticks = element_line(size = 0.8))  +
  theme(text = element_text(size = 13, color="black"))


#Plot
ggcorrplot(corr2018,type = "full", p.mat = p2018)+
  labs(x="",y="")+
  # theme_pubr()+
  # theme_article()+
  # mytheme+
  # theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text.x = element_text(size=14,angle = 45, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"))+
  theme(text = element_text(size = 16, color="black"))

