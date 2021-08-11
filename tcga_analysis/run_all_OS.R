remove(list = ls())

library(textreg)
library(tibble)
library(dplyr)
library(htmltab)
library(gridExtra)
library(mclust)
library(Cairo)
library(readr)
library(doParallel)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(ggplot2)
library(msigdbr)
library(tidyr)
library(gridExtra)
library(caret)
library(umap)
library(ggfortify)
library(mclust)
library(dbscan)
library(RColorBrewer)
library(arsenal)
library(dplR)
library(stargazer)
library(xtable)

setwd("./OS_Pipeline/R")

source("./OS_model_functions.R")

# get all the geneset combos
# relevant directories: DataRDS
dir.create("../DataRDS")
dir.create("../DataRDS/gene_sets/")
source("./saveSets.R")

# get 1D UMAP embeddings of all the patients
# relevant directories: UMAP_DataRDS, UMAP_1D
dir.create("../DataRDS/UMAP_embeddings/")
source("./UMAP_1D_embed.R")

dir.create("../DataTables")
dir.create("../DataTables/DBSCAN_Clust/") # the DBSCAN clustering for the KM model
dir.create("../DataTables/CoxPh/") # latex tables for Cox
dir.create("../DataTables/KM/") # latex tables for KM
dir.create("../Plots")
dir.create("../Plots/Surv_Plots/") # the plot of the DBSCAN clustering and corresponding KM model
source("./1D_Surv.R")
setwd("../../")
