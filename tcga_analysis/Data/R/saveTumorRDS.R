
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

for(i in tumortypes){ # gdc might throttle or something if we do full parallel here
  if(!file.exists(paste0("../raw_data/", i, "_data.RDS"))){
    data <- gdcQ(paste0("TCGA-",i),datadir="../GDCdata") 
    saveRDS(data, paste0("../raw_data/", i, "_data.RDS"))
  }
}
