cl2 <- makeCluster(14)
registerDoParallel(cl2)

allThree = osPfi = osDfi = c()
# no caveats on OS, DFI, PFI
allThree = c("BLCA", # Urothelial Bladder Carcinoma 
             "CESC", # Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma 
             "COAD", # Colon Adenocarcinoma
             "ESCA", # Esophageal Carcinoma
             "HNSC", # Head-Neck Squamous Cell Carcinoma
             "KIRP", # Cervical Kidney renal papillary cell carcinoma
             "LIHC", # Liver Hepatocellular Carcinoma 
             "LUAD", # Lung Adenocarcinoma
             "LUSC", # Lung Squamous Cell Carcinoma
             "OV", # Ovarian cancer
             "PAAD", # Pancreatic ductal adenocarcinoma
             "SARC", # Sarcoma
             "STAD", # Stomach Adenocarcinoma
             "UCEC") # Uterine Corpus Endometrial Carcinoma
# no caveats on OS and PFI
# osPfi = c("GBM",
#           "KIRC",
#           "MESO",
#           "SKCM",
#           "UVM")

tumortypes = c(allThree, osPfi, osDfi)

datadir = "../../Data/primary_tumors/"
genesetdir = "../DataRDS/gene_sets/"
outputdir = "../DataRDS/UMAP_embeddings/"

for(tumortype in tumortypes){
  if(!file.exists(paste0(outputdir, tumortype, "_UMAP_embeddings.RDS"))){
    data <- readRDS(paste0(datadir,tumortype,"_primary.RDS"))
    sets <- readRDS(paste0(genesetdir,tumortype,"_genesets.RDS"))
    embedding = embedMultiSet(sets = sets,
                              data = data$data,
                              tumortype = tumortype,
                              outputdir = "../DataRDS/UMAP_embeddings/",
                              n_components = 1)
  }
}

# for(tumortype in tumortypes){
#   source("./Cox_model_functions.R")
#   set.seed(1)
#   datadir = "./UMAP_DataRDS/"
#   data <- readRDS(paste0(datadir,tumortype,"_primary.RDS"))
#   sets <- readRDS(paste0(datadir,tumortype,"_genesets.RDS"))
#   embedding = embedMultiSet(sets = sets,
#                             data = data$data,
#                             tumortype = tumortype,
#                             outputdir = "./UMAP_1D/",
#                             n_components = 1)
# }

stopCluster(cl2)
