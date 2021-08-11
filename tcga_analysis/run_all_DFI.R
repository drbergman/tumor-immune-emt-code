remove(list = ls())

library(gridExtra)
library(mclust)
library(Cairo)
library(dplyr)
library(caret)
library(umap)
library(Rtsne)
library(kernlab)
library(R.matlab)
library(stringr)

setwd("./DFI_Pipeline/R")

source("./surv_DFI_Bayes_Single_OneDFI_AddRecGene.R")

dir.create("../DFI_Class_Plots/")