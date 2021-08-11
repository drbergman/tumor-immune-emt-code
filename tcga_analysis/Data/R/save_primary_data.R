
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

datadir = "../raw_data/"
outputdir = "../primary_tumors/"
codesdir = "../"
# get the sample type codes from TCGA
codes = htmltab("https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes", which = 2)
saveRDS(codes, paste0(codesdir, "TCGA_Sample_Codes.RDS"))

for(i in 1:length(tumortypes)){
  set.seed(1)
  tumortype = tumortypes[i]
  if(!file.exists(paste0(outputdir,tumortype,"_primary.RDS"))){
    data <- readRDS(paste0(datadir,tumortype,"_data.RDS"))
    data_primary = getSamples(data_list = data,codes=codes)
    saveRDS(data_primary, paste0(outputdir,tumortype,"_primary.RDS"))
  }
}