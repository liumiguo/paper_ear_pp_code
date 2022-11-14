#======================================================
#2018-2019 gene infromation 
#======================================================
#A Venn diagram was used to show genes in 2018 and 2019
#Fig. 2A

#Set path
# setwd(choose.dir())
# getwd()

#Load needed package 
library(tidyverse)
library(readxl)
library(ggvenn)

#Read data from xlsx file
list.files(".")
readxl::excel_sheets("./Data/2018-2019dataInformation.xlsx")

D2018<-read_xlsx("./Data/2018-2019dataInformation.xlsx",sheet = "2018")
D2019<-read_xlsx("./Data/2018-2019dataInformation.xlsx",sheet = "2019")

#Obtain genes in two years
Gene_both<-intersect(D2018$`Material ID`,D2019$`Material ID`)
Gene_both<-as.data.frame(Gene_both)


#Plot Venn diagram
Data_list<-list("2018"=D2018$`Material ID`,
             "2019"=D2019$`Material ID`)

ggvenn(Data_list,c("2018","2019"),
       fill_color = c("blue", "red"),
       fill_alpha = 0.6,
       set_name_size = 10,
       stroke_color = "grey",
       stroke_alpha=0.5,
       stroke_size = 0.5,
       text_size=7,
       text_color = "white")

