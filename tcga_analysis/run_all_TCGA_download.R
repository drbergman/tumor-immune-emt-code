remove(list = ls())

library(textreg)
library(tibble)
library(dplyr)
library(htmltab)
library(gridExtra)
library(mclust)
library(Cairo)
library(dplyr)
library(readr)


setwd("./Data/R")

source("./support_functions.R")

# get the GDC data by running:
# relevant directories: DataRDS
dir.create("../raw_data")
source("./saveTumorRDS.R")

# get only the primary tumor from those patients
# relevant directories: DataRDS
dir.create("../primary_tumors/")
source("./save_primary_data.R")

setwd("../../")