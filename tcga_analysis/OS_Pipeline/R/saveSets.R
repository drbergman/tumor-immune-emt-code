cl <- makeCluster(8)
registerDoParallel(cl)

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
foreach(i = tumortypes) %dopar% {
#for(i in tumortypes){
  if(!file.exists(paste0(genesetdir, i, "_genesets.RDS"))){
    data = readRDS(paste0(datadir, i, "_primary.RDS"))
    saveSets(tumortype = i,
             genes = data$data_genes,
             savedir = genesetdir)
  }
}

stopCluster(cl)
