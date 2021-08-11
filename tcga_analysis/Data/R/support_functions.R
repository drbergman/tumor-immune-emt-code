library(SummarizedExperiment)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(ggplot2)
library(msigdbr)
library(tidyr)
library(gridExtra)
library(caret)
library(Cairo)
library(umap)
library(dbscan)
library(mclust)

saveData <- function(tumortype = "PAAD",
                     datadir = "./DataRDS/"){
  data <- gdcQ(paste0("TCGA-",tumortype))
  saveRDS(data, paste0(datadir, tumortype, ".RDS"))
  return(data)
}

gdcQ <- function(proj = "TCGA-PAAD",datadir = "../GDCdata"){
  ## get the expression and survival data
  query <- GDCquery(project = proj,
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  GDCdownload(query, method = "api", files.per.chunk = 10, directory = datadir)
  data <- GDCprepare(query,directory = datadir)
  
  indDead = intersect(which(tolower(data$vital_status) == "dead"), which(data$days_to_death > 0))
  indAlive = Reduce(`intersect`, list(which(tolower(data$vital_status) == "alive"), which(is.na(data$days_to_death)), which(data$days_to_last_follow_up > 0)))  
  indAll = c(indDead, indAlive)
  
  data_matrix <- assay(data, 'normalized_count'); data_matrix = data_matrix[,indAll]
  data_genes <- rownames(data_matrix)
  time = c(data$days_to_death[indDead], data$days_to_last_follow_up[indAlive])
  status <- c(data$vital_status[indDead], data$vital_status[indAlive])
  bcr_patient_barcode <- data$bcr_patient_barcode
  
  bcr_patient_barcode = bcr_patient_barcode[indAll]
  
  for(i in 1:length(status)){
    if(tolower(status[i]) == 'alive'){
      status[i] <- 0 #these patients have been censored
    } else {
      status[i] <- 1
    }
  }
  
  return(list(data = data_matrix,
              data_genes = data_genes,
              status = status,
              time = time,
              bcr_patient_barcode = bcr_patient_barcode))
}

getSamples <- function(datadir="DataRDS",data_list = NULL,codes=NULL){
  # for a for a list of data, filter out everything except the primary tumor
  patients = unlist(lapply(strsplit(colnames(data_list$data), "-"), function(x){
    x[3]
  }))
  
  samples = unlist(lapply(strsplit(colnames(data_list$data), "-"), function(x){
    x[4]
  }))
  
  sam.type = substr(samples, start = 1, stop = 2)
  sam.vial = substr(samples, start = 3, stop = 3)
  
  # find all the patients that have primary tumors available
  prim_ind = which(sam.type == "01")
  
  pat.uni = unique(unlist(lapply(strsplit(colnames(data_list$data[,prim_ind]), "-"), function(x){
    x[3]
  })))
  ind = c()
  for(i in 1:length(pat.uni)){
    pat_all = which(patients == pat.uni[i])
    primary_all = which(sam.type == "01")
    primary_ind = intersect(pat_all, primary_all)
    if(length(primary_ind) > 1){
      primary_ind = primary_ind[1]
    }
    if(length(primary_ind) == 0){
      next
    }
    ind[i] = primary_ind    
  }
  cat("\n selecting ", length(ind), " of ", length(prim_ind), " primary samples \n")
  primary_data = data_list$data[,ind]
  primary_data_genes = rownames(primary_data)
  primary_status = data_list$status[ind]
  primary_time = data_list$time[ind]
  primary_bcr_patient_barcode = data_list$bcr_patient_barcode[ind]
  
  return(list(data = primary_data,
              data_genes = primary_data_genes,
              status = primary_status,
              time = primary_time,
              bcr_patient_barcode = primary_bcr_patient_barcode))
}

