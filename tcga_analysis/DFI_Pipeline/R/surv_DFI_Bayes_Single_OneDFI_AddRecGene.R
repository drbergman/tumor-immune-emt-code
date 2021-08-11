OSresults_path = "../../OS_Pipeline/DataTables/OS_test_results_filtered.csv"
UMAP_path = "../../OS_Pipeline/DataRDS/UMAP_embeddings/"
data_path = "../../Data/primary_tumors/"
gmplotdir = "../DFI_Class_Plots/"
savepath = "../DataTables/"

trainrate = .5

TCGA.ExtraEndPoints <- read.csv("../Liu_2018_Tables/TCGA-ExtraEndPoints.csv", row.names=1)

# get the mesenchymal cell proliferation gene sets
prolif_list = readRDS("../GO_genesets/GO_MESENCHYMAL_CELL_PROLIFERATION.RDS") # these are prolif genes
prolif_list = c(prolif_list,c("FZD2","FZD4","FZD6","FZD8","ROR1","ROR2","RYK","LRP6")) # add
bothHR_all = read.csv(OSresults_path, sep="", stringsAsFactors=FALSE)
bothHR_all = as_tibble(bothHR_all)

tumor_type = c()
prolif_set = c()
combo_row = c()
combo_cluster = c()
n_data_points = c()
train_rate = c()
pred_acc = c()
umap_plot = c()
fname_embedding = c()
fname_expression = c()
fname_expression_genes = c()
dirname_gpplot = c()
fname_gpml_vars = c()
negLogLik = c()
embed_type = c()

counter = 1

for(ttype in unique(bothHR_all$tumortype)){
  # 1) get all the data
  data_list <- readRDS(paste0(data_path, ttype, "_primary.RDS")) # save this from the filtered primary tumor samples
  colnames(data_list$data) = data_list$bcr_patient_barcode

  # 2) get the data with DFI
  ind=match(substr(data_list$bcr_patient_barcode,start=1,stop=12),as.character(TCGA.ExtraEndPoints$bcr_patient_barcode))[!is.na(match(substr(data_list$bcr_patient_barcode,start=1,stop=12),as.character(TCGA.ExtraEndPoints$bcr_patient_barcode)))] # this is the index in the Liu table
  PFI_incl_all = ind[which(as.character(TCGA.ExtraEndPoints$PFI.cr[ind])=="1"|
                           as.character(TCGA.ExtraEndPoints$PFI.cr[ind])=="2"|
                           as.character(TCGA.ExtraEndPoints$PFI.cr[ind])=="0")]
  DFI_incl_all = ind[which(as.character(TCGA.ExtraEndPoints$DFI.cr[ind])=="1"|
                           as.character(TCGA.ExtraEndPoints$DFI.cr[ind])=="2"|
                           as.character(TCGA.ExtraEndPoints$DFI.cr[ind])=="0")]
  PFI_incl_cen = ind[which(as.character(TCGA.ExtraEndPoints$PFI.cr[ind])=="1"|
                           as.character(TCGA.ExtraEndPoints$PFI.cr[ind])=="0")]
  DFI_incl_cen = ind[which(as.character(TCGA.ExtraEndPoints$DFI.cr[ind])=="1"|
                           as.character(TCGA.ExtraEndPoints$DFI.cr[ind])=="0")]
  PFI_uncen = ind[which(as.character(TCGA.ExtraEndPoints$PFI.cr[ind])=="1")]
  DFI_uncen = ind[which(as.character(TCGA.ExtraEndPoints$DFI.cr[ind])=="1")]
  
  #subsetDFI <- select(TCGA.ExtraEndPoints[DFI_incl_cen,],c(bcr_patient_barcode,DFI.time.cr))
  subsetDFI <- select(TCGA.ExtraEndPoints[DFI_incl_all,],c(bcr_patient_barcode,DFI.time.cr))
  subsetDFI$DFI.time.cr <- as.numeric(subsetDFI$DFI.time.cr)
  
  # 3) get the expression matrix indices of the Liu endpoints
  data.bar = substr(data_list$bcr_patient_barcode, start = 9, stop = 12)
  dfi.bar = substr(subsetDFI$bcr_patient_barcode, start = 9, stop = 12)
  ind_dfi = match(dfi.bar, data.bar)[which(!is.na(match(dfi.bar, data.bar)))] # the clustering ids (filtered by dfi)
  
  fit = Mclust(subsetDFI$DFI.time.cr, G=2, model="V")
  dfi= fit$classification
  dfi = (dfi-1)*2-1 # make low DFI -1 and high DFI +1 for GP classification

  # get the genes from the gene set that we want to model, the loop filters genes with low expression
  rows <- match(prolif_list, rownames(data_list$data))
  rows2 = c()
  for(k in 1:length(rows)){
    print(paste0(sum(data_list$data[rows[k],ind_dfi] == 0), " of ", length(data_list$data[rows[k],ind_dfi]), " are zero for ", 
                 rownames(data_list$data)[rows[k]], " of ", ttype))
    if(sum(data_list$data[rows[k],ind_dfi] == 0) < 0.5*length(data_list$data[rows[k],ind_dfi])){
      rows2 = c(rows2, rows[k]) # only include the genes that are expressed in at least 10% of the patient tumors.  this may improve stability for the GPML
    }
  }
  rows = rows2
  data_prolif <- data_list$data[rows,ind_dfi]
  classdata = as.data.frame(t(data_prolif))
  # save the non-embedded data and response for GPML
  fname_expr <- paste0(ttype, "_EXPR_AddRecGene.csv")
  write.table(data.frame(log(classdata+1), dfi, barcodes = dfi.bar), paste0(savepath,fname_expr), col.names = F, row.names = F)
  # matlab does not like to import row and column names, so save the genes separately (in general, these are not the same between different rows of bothHR_all ie different cluster/combo/tumortypes)
  fname_genelist <- paste0(ttype, "_genelist_AddRecGene.csv")
  write.table(data.frame(colnames = c(rownames(data_list$data)[rows], "dfi", "barcode"), # all the genes, the class and the barcode
                         datatypes = c(rep("double",length(rows)), "double", "char")), # the types (for readtable matlab options)  
              paste0(savepath,fname_genelist), row.names = F)
  
  fname <- paste0(gmplotdir, ttype, "_DFI_class.pdf")
  pdf(file = fname, width = 5, height = 8)
  par(mfrow=c(3,1))
  plot(fit, what = c("classification"),xlab = "Disease Free Interval (Days)")
  plot(fit, what = c("density"),xlab = "Disease Free Interval (Days)")
  hist(subsetDFI$DFI.time.cr,xlab = "Disease Free Interval (Days)", main = NULL)
  dev.off()
  
  tumor_type[counter] = ttype
  n_data_points[counter] = ncol(data_prolif)
  train_rate[counter] = trainrate
  pred_acc[counter] = 0
  fname_expression[counter] = fname_expr
  fname_expression_genes[counter] = fname_genelist
  dirname_gpplot[counter] = paste0(ttype)
  counter = counter + 1
}

Prolif_acc = data.frame(tumor_type,n_data_points,pred_acc, fname_expression, fname_expression_genes, dirname_gpplot) # each row of this table will be one GPML model fitted by running gp_tumor_single.m
write.table(Prolif_acc, "../DataTables/Prolif_acc_AddRecGene.txt", sep = "\t", row.names = F)
