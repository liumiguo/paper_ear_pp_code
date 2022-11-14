#========================================================================
#Distribution of the G,B,E value 2019
#========================================================================
# setwd("./Data/")
getwd()

source("./Code/Tool/00.load-packages.R")
library(ggpubr)

# Load and prep data ------------------------------------------------------
# traitMatrix <- readRDS("data/tidy_traitMatrix.rds")
earPhe<-read_rds("./Data/ear_trait_data_2019.rds")
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
  new.res <- readRDS(paste0("./Data/2019fwr-results/", phenos[i], "_IQR-fwr.rds"))
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

#calculating the distribution
# G value
names(g.table)

g.table2<-g.table%>%
  within(.,{
    Phenotype2<-NA
    Phenotype2[Phenotype=="MBaldTip"]<-"Bald tip length\nGenetic value"
    Phenotype2[Phenotype=="MBaldTipAreaRatio"]<-"Proportion of bald tip area\nGenetic value"
    Phenotype2[Phenotype=="MDiseaseAreaRatio"]<-"Proportion of disease area\nGenetic value"
    Phenotype2[Phenotype=="MEarCircumference"]<-"Ear circumference\nGenetic value"
    Phenotype2[Phenotype=="MEarDiameter"]<-"Ear diameter\nGenetic value"
    Phenotype2[Phenotype=="MEmptyAreaRatio"]<-"Proportion of empty area\nGenetic value"
    Phenotype2[Phenotype=="MNormalSeedAreaRatio"]<-"Proportion of normal kernel area\nGenetic value"
    Phenotype2[Phenotype=="MNormalSeedNum"]<-"Normal kernel number\nGenetic value"
    Phenotype2[Phenotype=="MRowSeedNum"]<-"Number of kernels per row\nGenetic value"
    Phenotype2[Phenotype=="MSeedThickness"]<-"Kernel thickness\nGenetic value"
    Phenotype2[Phenotype=="MSeedWith"]<-"Kernel width\nGenetic value"
    Phenotype2[Phenotype=="MShrinkageAreaRatio"]<-"Proportion of shrink area\nGenetic value"
    Phenotype2[Phenotype=="MShrunkenSeedNum"]<-"Shrunken kernel number\nGenetic value"
    Phenotype2[Phenotype=="MSpikeL"]<-"Ear length\nGenetic value"
    Phenotype2[Phenotype=="MSpikeRow"]<-"Number of rows per ear\nGenetic value"
    Phenotype2[Phenotype=="MSpikeWidth"]<-"Ear width\nGenetic value"
    Phenotype2[Phenotype=="MSpikeShape"]<-"Ear shape\nGenetic value"
  })%>%
  dplyr::rename(value=G)

# g.table_wt<-g.table2%>%
#   filter(Genotype=="9999901001")
# P1<-ggplot(data = g.table2,aes(x=G))+
#   geom_density(fill="red",color="white",alpha=0.6)+
#   geom_vline(data = g.table_wt,aes(xintercept = G),color="blue",linetype=2)+
#   # geom_histogram()+
#   facet_wrap(Phenotype2~.,scales="free",nrow = 3,ncol = 5)+
#   # scale_y_sqrt()+
#   labs(title=" ",x="value")+
#   theme_bw()+
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())


# B value
names(b.table)

b.table2<-b.table%>%
  within(.,{
    Phenotype2<-NA
    Phenotype2[Phenotype=="MBaldTip"]<-"Bald tip length\nLinear plasticity"
    Phenotype2[Phenotype=="MBaldTipAreaRatio"]<-"Proportion of bald tip area\nLinear plasticity"
    Phenotype2[Phenotype=="MDiseaseAreaRatio"]<-"Proportion of disease area\nLinear plasticity"
    Phenotype2[Phenotype=="MEarCircumference"]<-"Ear circumference\nLinear plasticity"
    Phenotype2[Phenotype=="MEarDiameter"]<-"Ear diameter\nLinear plasticity"
    Phenotype2[Phenotype=="MEmptyAreaRatio"]<-"Proportion of empty area\nLinear plasticity"
    Phenotype2[Phenotype=="MNormalSeedAreaRatio"]<-"Proportion of normal kernel area\nLinear plasticity"
    Phenotype2[Phenotype=="MNormalSeedNum"]<-"Normal kernel number\nLinear plasticity"
    Phenotype2[Phenotype=="MRowSeedNum"]<-"Number of kernels per row\nLinear plasticity"
    Phenotype2[Phenotype=="MSeedThickness"]<-"Kernel thickness\nLinear plasticity"
    Phenotype2[Phenotype=="MSeedWith"]<-"Kernel width\nLinear plasticity"
    Phenotype2[Phenotype=="MShrinkageAreaRatio"]<-"Proportion of shrink area\nLinear plasticity"
    Phenotype2[Phenotype=="MShrunkenSeedNum"]<-"Shrunken kernel number\nLinear plasticity"
    Phenotype2[Phenotype=="MSpikeL"]<-"Ear length\nLinear plasticity"
    Phenotype2[Phenotype=="MSpikeRow"]<-"Number of rows per ear\nLinear plasticity"
    Phenotype2[Phenotype=="MSpikeWidth"]<-"Ear width\nLinear plasticity"
    Phenotype2[Phenotype=="MSpikeShape"]<-"Ear shape\nLinear plasticity"
  })%>%
  dplyr::rename(value=B)


# b.table_wt<-b.table2%>%
#   filter(Genotype=="9999901001")

gb.table<-rbind(g.table2,b.table2)%>%
  within(.,{
    Phenotype2<-factor(Phenotype2,ordered = T,levels = c(unique(g.table2$Phenotype2),unique(b.table2$Phenotype2)))
  })

gb.table_wt<-gb.table%>%
  filter(Genotype=="9999901001")

#Fig. S8
ggplot(data = gb.table,aes(x=value))+
  geom_density(fill="red",color="white",alpha=0.6)+
  geom_vline(data = gb.table_wt,aes(xintercept = value),color="blue",linetype=2)+
  # geom_histogram()+
  facet_wrap(Phenotype2~.,scales="free",nrow = 6,ncol = 5)+
  # scale_y_sqrt()+
  labs(title=" ",x="Value",y="Density")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())


names(gb.table)
names(g.table2)
gb.table2<-left_join(g.table,b.table)%>%
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
names(gb.table2)

gb.table2_wt<-gb.table2%>%
  filter(Genotype=="9999901001")

#Fig. S10
ggplot(data = gb.table2,aes(x=B,y=G))+
  geom_point(color="red",alpha=0.5)+
  # geom_histogram()+
  geom_vline(data=gb.table2_wt,aes(xintercept = B),color="blue")+
  # geom_vline(data=gb.table2_wt,aes(xintercept = 1),color="blue",linetype=2)+
  geom_hline(data=gb.table2_wt,aes(yintercept = G),color="blue")+
  facet_wrap(Phenotype2~.,scales="free",nrow = 6,ncol = 5)+
  # scale_y_sqrt()+
  labs(title=" ",x="Linear plasticity",y="Genetic value")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())


#---gbe相关性分析 Fig. S7

gbe.table<-Reduce(left_join,list(g.table,b.table,e.table))%>%
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
  })%>%
  select(Phenotype2,Genotype,G,B,VE)%>%
  rename(`Mean Phenotype`=G,
         `Linear Plasticity`=B,
         `Non-linear Plasticity`=VE)

phenos <- unique(gbe.table$Phenotype2)
traitList <- split(gbe.table, gbe.table$Phenotype2)


grobs1 <- 1:9 %>% map(function(i) {
  cordata<-traitList[[i]]%>%drop_na()
  corvalue<-cor(cordata[,c("Mean Phenotype","Linear Plasticity","Non-linear Plasticity")])

  test2 <- tibble(x = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), each = 3),
                  y = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), times = 3),
                  gcor = c(corvalue[, 1], corvalue[, 2], corvalue[, 3])) %>%
    mutate(x = ordered(x, levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity")),
           y = ordered(y, levels = c("Non-linear Plasticity", "Linear Plasticity", "Mean Phenotype")))
  ggplot(test2, aes(x = x, y = y, fill = gcor)) +
    geom_tile() + theme_bw() + labs(x = "", y = "", fill = expression(r[g])) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(phenos[i])

})

library(gridExtra)
grid.arrange(grobs = grobs1, ncol = 3)


grobs2 <- 10:15 %>% map(function(i) {
  cordata<-traitList[[i]]%>%drop_na()
  corvalue<-cor(cordata[,c("Mean Phenotype","Linear Plasticity","Non-linear Plasticity")])

  test2 <- tibble(x = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), each = 3),
                  y = rep(c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity"), times = 3),
                  gcor = c(corvalue[, 1], corvalue[, 2], corvalue[, 3])) %>%
    mutate(x = ordered(x, levels = c("Mean Phenotype", "Linear Plasticity", "Non-linear Plasticity")),
           y = ordered(y, levels = c("Non-linear Plasticity", "Linear Plasticity", "Mean Phenotype")))
  ggplot(test2, aes(x = x, y = y, fill = gcor)) +
    geom_tile() + theme_bw() + labs(x = "", y = "", fill = expression(r[g])) +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(phenos[i])

})

library(gridExtra)
grid.arrange(grobs = grobs2, ncol = 3)






# P2<-ggplot(data = b.table2,aes(x=B))+
#   geom_density(fill="red",color="white",alpha=0.6)+
#   geom_vline(data = b.table_wt,aes(xintercept = B),color="blue",linetype=2)+
#   # geom_histogram()+
#   facet_wrap(Phenotype2~.,scales="free",nrow = 3,ncol = 5)+
#   # scale_y_sqrt()+
#   labs(title=" ",x="value")+
#   theme_bw()+
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())


# ggpubr::ggarrange(P1,P2,nrow = 2)


# e value
# names(e.table)
# 
# e.table2<-e.table%>%
#   within(.,{
#     Phenotype2<-NA
#     Phenotype2[Phenotype=="MBaldTip"]<-"Bald tip length\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MBaldTipAreaRatio"]<-"The proportion of bald tip area\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MDiseaseAreaRatio"]<-"The proportion of disease area\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MEarCircumference"]<-"Ear circumference\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MEarDiameter"]<-"Ear diameter\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MEmptyAreaRatio"]<-"The proportion of empty area\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MNormalSeedAreaRatio"]<-"The proportion of normal kernels area\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MNormalSeedNum"]<-"Normal kernels number\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MRowSeedNum"]<-"The number of kernels per row\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MSeedThickness"]<-"Kernels thickness\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MSeedWith"]<-"Kernels Width\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MShrinkageAreaRatio"]<-"The proportion of shrink area\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MShrunkenSeedNum"]<-"Shrunken kernels number\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MSpikeL"]<-"Ear length\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MSpikeRow"]<-"The number of rows per ear\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MSpikeWidth"]<-"Ear width\nNon-Linear plasticity"
#     Phenotype2[Phenotype=="MSpikeShape"]<-"Ear shape\nNon-Linear plasticity"
#   })
# 
# e.table_wt<-e.table2%>%
#   filter(Genotype=="9999901001")
# 
# ggplot(data = e.table2,aes(x=VE))+
#   geom_density(fill="red",color="white",alpha=0.6)+
#   geom_vline(data = e.table_wt,aes(xintercept = VE),color="blue",linetype=2)+
#   # geom_histogram()+
#   facet_wrap(Phenotype2~.,scales="free",nrow = 3,ncol = 5)+
#   # scale_y_sqrt()+
#   labs(title="2018",x="value")+
#   theme_bw()+
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank())
