#==========================================================
#maize ear phenotypes of inbred line in 2018
#==========================================================

# setwd("D:/Post_Doctor_Data_Analysis/Data/(整理)ND101-2014-2020田测数据/整理后的数据")
# getwd()
### Run Bayesian Finlay-Wilkinson regression for all 15 phenotype.
### This code assumes that you have 23 cores to run all regressions in parallel.
### my CP have 9 cores

datapath<-"./Data/ear_trait_data_2018.rds"#data file path

source("./Code/Tool/00.load-packages.R")
# library(devtools)
# install_github("lian0090/FW")
library(parallel)
library(doParallel)
library(foreach)
library(iterators)

# Set-up control variables for parallel processing
ncores <- detectCores()
cl <- makeCluster(ncores, outfile = "./Data/2018logs/bayesian_fwr_IQR.log")
registerDoParallel(cl)

# Load the trait matrix and split into a list
# traitMatrix <- readRDS("data/tidy_traitMatrix_IQR_AK.rds")
earPhe<-read_rds(datapath)
names(earPhe)
earPhe<-earPhe%>%
  mutate(Environment=str_c(Site,"_",Rep))

# phenos <- unique(earPhe$Index)
traitList <- split(earPhe, earPhe$Index)
phenos<-names(traitList)
### Bayesian FWR in parallel
# Control variables
burnin <- 1000
niter <- 51000

# Random seeds
seeds <- 47060859 + seq(0, ncores - 1)

fw <- foreach(ph = 1:15, .packages = 'FW') %dopar% {
  result <- with(traitList[[ph]],
                 FW(y = value, VAR = Genotype, ENV = Environment, seed = seeds[ph],
                    method = "Gibbs", nIter = niter, burnIn = burnin,
                    saveAt = paste0("./Data/2018gibbs-samples/", phenos[ph], "_IQR-Gibbs")))
  saveRDS(result, paste0("./Data/2018fwr-results/", phenos[ph], "_IQR-fwr.rds"))
}

stopCluster(cl)
