library(tidyverse)
library(stringr)
library(broom)

# Colors for plotting results
colors <- c("royalblue", "orange", "red")

# Vector of phenotype names
phenos <- c("MBaldTip", "MEarCircumference", "MEarDiameter", "MSpikeL", "MSpikeShape","MSeedThickness",
            "MSpikeWidth", "MSeedWith", "MNormalSeedNum", "MShrunkenSeedNum", "MRowSeedNum",
            "MSpikeRow", "MBaldTipAreaRatio", "MEmptyAreaRatio",
            "MNormalSeedAreaRatio")

# Vector pretty phenotype names for plotting
lbls <- c("Bald tip length (mm)","Ear circumference (mm)","Ear diameter (mm)","Ear length (mm)",
          "Ear shape","Ear width (mm)","Kernels thickness (mm)",
          "Kernels Width (mm)","Normal seed number","Shrunken seed number",
          "The number of kernels per row","The number of row per ear",
          "The proportion of bald tip area","The proportion of empty area","The proportion of normal seed area")
