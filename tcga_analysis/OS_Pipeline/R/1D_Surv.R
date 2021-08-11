source("./OS_model_functions.R")
saveall=FALSE # set this to true to save all the plots and embeddings regardless of Coxph test

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
# # no caveats on OS and PFI
# osPfi = c("GBM",
#           "KIRC",
#           "MESO",
#           "SKCM",
#           "UVM")

tumortypes = c(allThree, osPfi, osDfi)

datadir = "../../Data/primary_tumors/"
genesetdir = "../DataRDS/gene_sets/"
umapdir = "../DataRDS/UMAP_embeddings/"

alltumors = tibble()
for(tumortype in tumortypes){
  set.seed(1)
  
  data <- readRDS(paste0(datadir,tumortype,"_primary.RDS"))

  sets <- readRDS(paste0(genesetdir,tumortype,"_genesets.RDS"))

  embedding = readRDS(paste0(umapdir, tumortype, "_UMAP_embeddings.RDS"))

  pv = get1Dpvals(tumortype = tumortype,data = data,sets = sets, embedding = embedding)
  pvtib = as_tibble(pv$summary)

  pvtib$Schoenfeld = unlist(lapply(pv$models, function(x){
    check = cox.zph(x)
    length(which(check$table[,"p"] <= 5e-2))
  })) # 0 indicates non-significant relationship between Schoenfeld residuals and time
  pv$summary$Schoenfeld = unlist(lapply(pv$models, function(x){
    check = cox.zph(x)
    length(which(check$table[,"p"] <= 5e-2))
  })) # 0 indicates non-significant relationship between Schoenfeld residuals and time
  
  coxph_filter = c()
  if(dim(pvtib)[1] > 0){
    pvtib$tumortype = rep(tumortype, dim(pvtib)[1])
    pchsq_km = rep(Inf, dim(pvtib)[1])
    
    for(i in 1:nrow(pvtib)){
      coxph_filter[i] = (pv$summary$ScoreTest[i] <= 5e-2 &
                        pv$summary$Schoenfeld[i] == 0 &
                        pv$summary$bothP[i] <=5e-2 & pv$summary$emtP[i] <= 5e-2 & pv$summary$inflamP[i] <= 5e-2 &
                        pv$summary$inflamHR[i] < pv$summary$bothHR[i] & pv$summary$emtHR[i] < pv$summary$bothHR[i] &
                        pv$summary$bothHR[i] >= (pv$summary$emtHR + 0.05)[i] &
                        pv$summary$bothHR[i] >= (pv$summary$inflamHR[i] + 0.05))
      
      if(coxph_filter[i]+saveall){
        print(paste0("checking ", tumortype, " pvtib, row ", i, " of ", nrow(pvtib)))
        pchsq_km[i] = plotTT(embedding = embedding$BOTH$embed,
                                   models = pv,
                                   tumortype = tumortype,
                                   plotdir = "../Plots/Surv_Plots/",
                                   comboInd = i,
                                   data = data)
      }      
    }
  }

  pvtib$pchsq_km = pchsq_km
  pvtib$coxph_filter = coxph_filter
  alltumors = bind_rows(alltumors,pvtib)
}
browser()
write.table(alltumors, paste0("../DataTables/OS_test_results.csv"), row.names = F)

alltumors_filtered = filter(alltumors,pchsq_km <= 5e-2) # filtered table for later use
write.table(alltumors_filtered, paste0("../DataTables/OS_test_results_filtered.csv"), row.names = F)




