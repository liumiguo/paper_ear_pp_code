#========================================================================
#quartile coefficient of dispersion of linear and non-linear plasticity
#========================================================================
# setwd()
getwd()
source("./Code/Tool/00.load-packages.R")
library(data.table)

# Load and prep data ------------------------------------------------------
# traitMatrix <- readRDS("data/tidy_traitMatrix.rds")
earPhe<-read_rds("./Data/ear_trait_data_2018.rds")
names(earPhe)
earPhe<-earPhe%>%
  mutate(Environment=str_c(Site,"_",Rep))

phenos <- unique(earPhe$Index)
rm(earPhe)

# Create tables to hold results
g.table <- tibble(Phenotype = "", Genotype = "", G = 0, Set = "")
b.table <- tibble(Phenotype = "", Genotype = "", B = 0, Set = "")
e.table <- tibble(Phenotype = "", Genotype = "", VE = 0, Set = "")


# Main loop ---------------------------------------------------------------
# Main loop to gather calculated values
for (i in seq_along(phenos)) {
  # Load the FW results
  # i<-1
  new.res <- readRDS(paste0("./Data/2018fwr-results/", phenos[i], "_IQR-fwr.rds"))
  # load(paste0("2018gibbs-samples/", phenos[i], "_IQR-Gibbssamps.rda"))
  # load(paste0("~/gxe-gwas/data/", phenos[i], "2-gibbs.RData"))

  # Calculate residual variances
  var.new <- tibble(Genotype = new.res$VAR, y = new.res$y, yhat = new.res$yhat[, 1]) %>%
    mutate(e = y - yhat) %>%
    group_by(Genotype) %>%
    dplyr::summarise(VE = var(e)) %>%
    mutate(Set = "New", Phenotype = phenos[i]) %>%
    select(Phenotype, Genotype, VE, Set)
  # var.old <- tibble(Genotype = samps$VAR, y = samps$y, yhat = samps$yhat[, 1]) %>%
  #   mutate(e = y -yhat) %>%
  #   group_by(Genotype) %>%
  #   dplyr::summarise(VE = var(e)) %>%
  #   mutate(Set = "Old", Phenotype = phenos[i]) %>%
  #   select(Phenotype, Genotype, VE, Set)

  g.table <- bind_rows(list(g.table,
                            tibble(Phenotype = phenos[i], Genotype = new.res$VARlevels,
                                   G = new.res$g[, 1], Set = "New")))
  b.table <- bind_rows(list(b.table,
                            tibble(Phenotype = phenos[i], Genotype = new.res$VARlevels,
                                   B = new.res$b[, 1] + 1, Set = "New")))
  e.table <- bind_rows(list(e.table, var.new))
}

# Remove the initial value in each tibble
g.table <- g.table[-1, ]
b.table <- b.table[-1, ]
e.table <- e.table[-1, ]


names(b.table)

b.table_summ<-b.table%>%
  group_by(Phenotype)%>%
  summarise(QCD_B=(quantile(B, probs = 0.75)-quantile(B, probs = 0.25))/(quantile(B, probs = 0.75)+quantile(B, probs = 0.25)))

names(e.table)
e.table_summ<-e.table%>%
  drop_na()%>%
  group_by(Phenotype)%>%
  summarise(QCD_E=(quantile(VE, probs = 0.75)-quantile(VE, probs = 0.25))/(quantile(VE, probs = 0.75)+quantile(VE, probs = 0.25)))

QCD_DATA<-full_join(b.table_summ,e.table_summ)%>%
  dplyr::rename(`Linear plasticity`=QCD_B,
                `Non-linear plasticity`=QCD_E)%>%
  gather(.,key = Index,value = value,`Linear plasticity`:`Non-linear plasticity`)%>%
  within(.,{
    Phenotype2<-NA
    Phenotype2[Phenotype=="MBaldTip"]<-"Bald tip length"
    Phenotype2[Phenotype=="MBaldTipAreaRatio"]<-"Proportion of bald tip area"
    Phenotype2[Phenotype=="MDiseaseAreaRatio"]<-"Proportion of disease area"
    Phenotype2[Phenotype=="MEarCircumference"]<-"Ear circumference"
    Phenotype2[Phenotype=="MEarDiameter"]<-"Ear diameter"
    Phenotype2[Phenotype=="MEmptyAreaRatio"]<-"Proportion of empty area"
    Phenotype2[Phenotype=="MNormalSeedAreaRatio"]<-"Proportion of normal kernel area"
    Phenotype2[Phenotype=="MNormalSeedNum"]<-"Normal kernel number"
    Phenotype2[Phenotype=="MRowSeedNum"]<-"Number of kernels per row"
    Phenotype2[Phenotype=="MSeedThickness"]<-"Kernel thickness"
    Phenotype2[Phenotype=="MSeedWith"]<-"Kernel width"
    Phenotype2[Phenotype=="MShrinkageAreaRatio"]<-"Proportion of shrink area"
    Phenotype2[Phenotype=="MShrunkenSeedNum"]<-"Shrunken kernel number"
    Phenotype2[Phenotype=="MSpikeL"]<-"Ear length"
    Phenotype2[Phenotype=="MSpikeRow"]<-"Number of rows per ear"
    Phenotype2[Phenotype=="MSpikeWidth"]<-"Ear width"
    Phenotype2[Phenotype=="MSpikeShape"]<-"Ear shape"
  })


p2018<-ggplot(data = QCD_DATA,aes(x=Phenotype2,y=value,shape=Index,color=Index))+
  geom_point(size=5,alpha=0.7)+
  scale_color_manual(values = c("red","blue"))+
  labs(x=" ", y=" Quartile coefficient of dispersion\n",color=" ",shape=" ")+
  ylim(0,0.4)+
  # theme_pubr()+
  # theme_article()+
  theme_bw()+
  theme(legend.position = "top",
        axis.text.x = element_text(size=14,angle = 45, hjust = 1,color="black"),
        axis.text.y = element_text(size=14,color="black"),
        text = element_text(size=14,color="black"))

saveRDS(p2018,"./Data/2018_QCD2.rds")


names(g.table)
g.table_summ<-g.table%>%
  drop_na()%>%
  group_by(Phenotype)%>%
  summarise(QCD_E=(quantile(G, probs = 0.75)-quantile(G, probs = 0.25))/(quantile(G, probs = 0.75)+quantile(G, probs = 0.25)))

