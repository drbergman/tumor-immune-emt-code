## outputs a list of tables, one for each keyword set
## each row contains a set name and gene list
kwsearch <- function(data_genes,
                     species = 'Homo sapiens',
                     keywords =
                       list(EMT = c("EPITHELIAL_TO_MESENCHYMAL",
                                    "EMT"),
                            INFLAM = c("INFLAMMATORY",
                                       "INFLAMMATION"))){
  
  ## look for some interesting gene sets (axis 1 and 2)
  # load all msigdb sets
  #browser()
  msigdb <- msigdbr::msigdbr(species = species)
  set_list <- list()
  for(i in 1:length(keywords)){
    # find interesting sets
    key <- keywords[[i]]
    set_names <- unique(grep(paste(key,collapse="|"), 
                             msigdb$gs_name, value=TRUE))
    
    gs_table <- data.frame(matrix(0, length(set_names),0))
    
    ## add the genes to tables (axis 1 and 2)
    genes <- list()
    names <- c()
    for(j in 1:length(set_names)){
      # make sure we have data on the genes in this set
      ind <- grep(set_names[j], msigdb$gs_name)
      filter <- which(msigdb$gene_symbol[ind] %in% data_genes)
      ind <- ind[filter]
      names[j] <- set_names[j]
      genes[[j]] <- msigdb$gene_symbol[ind]
    }
    tbl_ind <- match(names, msigdb$gs_name)
    set_ids <- msigdb$gs_id[tbl_ind]
    
    gs_table$set_names <- names
    gs_table$set_ids <- set_ids
    gs_table$genes <- genes
    set_list[[as.name(names(keywords)[i])]] <- gs_table
  }
  return(set_list)
}

saveSets <- function(tumortype = "PAAD",
                     genes = NULL,
                     savedir = "./DataRDS/"){
  sets = list()
  print("searching BOTH")
  # save the gene sets for various keyword combos
  kw1 <- list(EMT = c("EPITHELIAL_TO_MESENCHYMAL",
                      "EMT"),
              INFLAM = c("INFLAMMATORY",
                         "INFLAMMATION"))
  sets[["both"]] <- kwsearch(keywords = kw1,
                             data_genes = genes)
  
  print("searching EMT")
  kw2 <- list(EMT = c("EPITHELIAL_TO_MESENCHYMAL","EMT"))
  sets[["emt"]] <- kwsearch(keywords = kw2,
                            data_genes = genes)
  
  
  print("searching INFLAM")
  kw3 <- list(INFLAM = c("INFLAMMATORY",
                         "INFLAMMATION"))
  sets[["inflam"]] <- kwsearch(keywords = kw3,
                               data_genes = genes)
  saveRDS(sets, paste0(savedir, tumortype, "_genesets.RDS"))
}

get1Dpvals <- function(tumortype = NULL,
                     data_list = NULL,
                     sets = NULL,
                     embedding = embedding,
                     cols = NULL){
  # get the log rank test pvals for each gene list combo of the given subset of the data
  pvals = list()
  
  if(is.null(cols)){
    cols = 1:ncol(data_list$data) # if this is not a fold then use all the data
  }
  
  models = list()
  for(i in 1:length(embedding$BOTH$embed)){
    emt_ind = which(embedding$EMT$names == embedding$BOTH$names_EMT[i])
    inflam_ind = which(embedding$INFLAM$names == embedding$BOTH$names_INFLAM[i])
    dataMerge = as.data.frame(cbind(embedding$EMT$embed[[emt_ind]],
                                    embedding$INFLAM$embed[[inflam_ind]],
                                    embedding$BOTH$embed[[i]]))
    colnames(dataMerge) = c("EMT", "INFLAM", "BOTH")
    dataMerge$status = as.numeric(data_list$status)
    dataMerge$time = as.numeric(data_list$time)
    f = as.formula(
      paste("Surv(time, status)", 
            paste(c("EMT", "INFLAM", "BOTH"), 
                  collapse = " + "), 
            sep = " ~ "))
    
    models[[i]] = coxph(f, data = dataMerge)
  }
  
  hr_log_summary = data.frame(matrix(rep(0,(11*length(models))),ncol = 11))
  colnames(hr_log_summary) = c("Row", "SetEMT", "SetINFLAM",
                               "emtHR", 
                               "inflamHR", 
                               "bothHR",
                               "emtP", 
                               "inflamP", 
                               "bothP",
                               "LikTest",
                               "ScoreTest")
  for(i in 1:length(models)){
    hr_log_summary$Row[i] = i
    hr_log_summary$SetEMT[i] = embedding$BOTH$names_EMT[i]
    hr_log_summary$SetINFLAM[i] = embedding$BOTH$names_INFLAM[i]
    
    hr_log_summary$emtHR[i] = summary(models[[i]])$coefficients[1,2]
    hr_log_summary$inflamHR[i] = summary(models[[i]])$coefficients[2,2]
    hr_log_summary$bothHR[i] = summary(models[[i]])$coefficients[3,2]

    hr_log_summary$emtP[i] = summary(models[[i]])$coefficients[1,5]
    hr_log_summary$inflamP[i] = summary(models[[i]])$coefficients[2,5]
    hr_log_summary$bothP[i] = summary(models[[i]])$coefficients[3,5]

    hr_log_summary$LikTest[i] = summary(models[[i]])$logtest[3]
    hr_log_summary$ScoreTest[i] = summary(models[[i]])$sctest[3]
  }
  
  return(list(summary = hr_log_summary, models = models))
}

getPvals <- function(tumortype = NULL,
                     data_list = NULL,
                     sets = NULL,
                     embedding = embedding,
                     cols = NULL){
  # get the log rank test pvals for each gene list combo of the given subset of the data
  pvals = list()

  if(is.null(cols)){
    cols = 1:ncol(data_list$data) # if this is not a fold then use all the data
  }
  
  models = list()
  for(i in 1:length(embedding$BOTH$embed)){
    emt_ind = which(embedding$EMT$names == embedding$BOTH$names_EMT[i])
    inflam_ind = which(embedding$INFLAM$names == embedding$BOTH$names_INFLAM[i])
    dataMerge = as.data.frame(cbind(embedding$EMT$embed[[emt_ind]],
                                    embedding$INFLAM$embed[[inflam_ind]],
                                    embedding$BOTH$embed[[i]]))
    colnames(dataMerge) = c("EMT_01", "EMT_02", "INFLAM_01", "INFLAM_02", "BOTH_01", "BOTH_02")
    dataMerge$status = as.numeric(data_list$status)
    dataMerge$time = as.numeric(data_list$time)
    f = as.formula(
      paste("Surv(time, status)", 
            paste(c("EMT_01", "EMT_02", 
                    "INFLAM_01", "INFLAM_02", 
                    "BOTH_01", "BOTH_02"), 
                  collapse = " + "), 
            sep = " ~ "))
    
    models[[i]] = coxph(f, data = dataMerge)
  }

  hr_log_summary = data.frame(matrix(rep(0,(21*length(models))),ncol = 21))
  colnames(hr_log_summary) = c("Row", "SetEMT", "SetINFLAM",
                               "emtHR01", "emtHR02", 
                               "inflamHR01", "inflamHR02", 
                               "bothHR01", "bothHR02",
                               "emtP01", "emtP02", 
                               "inflamP01", "inflamP02", 
                               "bothP01", "bothP02",
                               "LikTest", "SigLikMaxHR",
                               "ScoreTest","SigScoreMaxHR",
                               "AllSig", "AvgRatio")
  hr_log_summary$SigLikMaxHR = as.logical(hr_log_summary$SigLikMaxHR)
  hr_log_summary$SigScoreMaxHR = as.logical(hr_log_summary$SigScoreMaxHR)
  hr_log_summary$AllSig = as.logical(hr_log_summary$AllSig)
  for(i in 1:length(models)){
    EMTMaxInd = which.max(abs(summary(models[[i]])$coefficients[1:2,2]))
    INFLAMMaxInd = which.max(abs(summary(models[[i]])$coefficients[3:4,2]))
    bothMaxInd = which.max(abs(summary(models[[i]])$coefficients[5:6,2]))
    
    emtMax = summary(models[[i]])$coefficients[EMTMaxInd,2]
    inflamMax = summary(models[[i]])$coefficients[INFLAMMaxInd,2]
    bothMax = summary(models[[i]])$coefficients[bothMaxInd,2]
    
    emtMaxP = summary(models[[i]])$coefficients[EMTMaxInd,5]
    inflamMaxP = summary(models[[i]])$coefficients[INFLAMMaxInd,5]
    bothMaxP = summary(models[[i]])$coefficients[bothMaxInd,5]
    
    AllSig = emtMaxP<=5e-2 & inflamMaxP<=5e-2 & bothMaxP<=5e-2
    
    # get the avg absolute value of all EMT or INFLAM HRs that are significant
    ind_sum = 0
    n = 0
    for(j in 1:4){
      if(summary(models[[i]])$coefficients[j,5] <= 5e-2){
        ind_sum = ind_sum + abs(summary(models[[i]])$coefficients[j,2])
        n = n + 1
      }
    }
    if(n > 0){
      avg_ind = ind_sum/n
    }else{
      avg_ind = 0
    }

    # get the avg absolute value of all BOTH HRs that are significant
    both_sum = 0
    n = 0
    for(j in 5:6){
      if(summary(models[[i]])$coefficients[j,5] <= 5e-2){
        ind_sum = ind_sum + abs(summary(models[[i]])$coefficients[j,2])
        n = n + 1
      }
    }
    if(n > 0){
      avg_both = both_sum/n
    }else{
      avg_both = 0
    }
    
    if(avg_ind > 0 & avg_both > 0){
      avg_ratio = avg_both/avg_ind
    }else{
      avg_ratio = 0
    }
    
    bothMaxHR = 3 == which.max(c(emtMax, inflamMax, bothMax))
    logTestSig = summary(models[[i]])$logtest[3] <= 5e-2
    scoreTestSig = summary(models[[i]])$sctest[3] <= 5e-2
    SigLikMaxHR = bothMaxHR & logTestSig & bothMaxP<=5e-2
    SigScoreMaxHR = bothMaxHR & scoreTestSig & bothMaxP<=5e-2
    
    hr_log_summary$Row[i] = i
    hr_log_summary$SetEMT[i] = embedding$BOTH$names_EMT[i]
    hr_log_summary$SetINFLAM[i] = embedding$BOTH$names_INFLAM[i]
    
    hr_log_summary$emtHR01[i] = summary(models[[i]])$coefficients[1,2]
    hr_log_summary$emtHR02[i] = summary(models[[i]])$coefficients[2,2]
    hr_log_summary$inflamHR01[i] = summary(models[[i]])$coefficients[3,2]
    hr_log_summary$inflamHR02[i] = summary(models[[i]])$coefficients[4,2]
    hr_log_summary$bothHR01[i] = summary(models[[i]])$coefficients[5,2]
    hr_log_summary$bothHR02[i] = summary(models[[i]])$coefficients[6,2]
    
    hr_log_summary$emtP01[i] = summary(models[[i]])$coefficients[1,5]
    hr_log_summary$emtP02[i] = summary(models[[i]])$coefficients[2,5]
    hr_log_summary$inflamP01[i] = summary(models[[i]])$coefficients[3,5]
    hr_log_summary$inflamP02[i] = summary(models[[i]])$coefficients[4,5]
    hr_log_summary$bothP01[i] = summary(models[[i]])$coefficients[5,5]
    hr_log_summary$bothP02[i] = summary(models[[i]])$coefficients[6,5]
    
    hr_log_summary$LikTest[i] = summary(models[[i]])$logtest[3]
    hr_log_summary$SigLikMaxHR[i] = SigLikMaxHR
    hr_log_summary$ScoreTest[i] = summary(models[[i]])$sctest[3]
    hr_log_summary$SigScoreMaxHR[i] = SigScoreMaxHR
    hr_log_summary$AllSig[i] = AllSig
    hr_log_summary$AvgRatio[i] = avg_ratio
  }

  return(list(summary = hr_log_summary, models = models))
}

plotTT <- function(data_list = NULL,
                   tumortype = NULL,
                   embedding = NULL,
                   models = NULL,
                   plotdir = "../Plots/",
                   kmtabdir = "../DataTables/KM/",
                   coxtabdir = "../DataTables/CoxPh/",
                   comboInd = NULL,
                   mingroupsize = 30,
                   dbscan_dir = "../DataTables/DBSCAN_Clust/"){
  # Save model summaries for all coxph 
  # If there are significant KM models, save the plots and the model summaries
  
  
  #######################
  ### save coxph latex table
  #######################
  captext = paste0("Summary of Cox PH Model for ", tumortype)
  fname = paste0(coxtabdir, tumortype, "_", comboInd, "_cph", ".tex")
  cat('\n saving ',fname,'\n')
  foo = capture.output(stargazer(models$models[[comboInd]], 
                                 label = paste0("table:", tumortype, "_cph"),
                                 title = latexify(captext, doublebackslash = F), # can set this to captext but seems to be creating separate tables
                                 out = fname))
  plots = list()
  colors = c("#b2df8a",
             "#33a02c")
  
  title="A. DBSCAN Clustering"
  set.seed(1)
  dbclust = hdbscan(x = embedding[[comboInd]], minPts = mingroupsize)
  
  while(length(unique(dbclust$cluster)) < 2){
    mingroupsize = 0.90*mingroupsize
    set.seed(1)
    dbclust = hdbscan(x = embedding[[comboInd]], minPts = mingroupsize)
  }
  ncolors = length(which(unique(dbclust$cluster) != 0))
  colors = brewer.pal(ncolors, name="Greens")
  
  ind = which(dbclust$cluster > 0) # 0 indicates noise points
  df = as.data.frame(cbind(embedding[[comboInd]][ind,], dbclust$cluster[ind]))
  if(dim(as.matrix(embedding[[comboInd]][ind,]))[2] == 2){
    colnames(df) = c("U1", "U2", "Cluster")
    df$Cluster = as.factor(df$Cluster)
    plots$scat.plot <-
      ggplot(df,
             aes(x = U1,
                 y = U2,
                 color = Cluster)) +
      geom_point() +
      scale_color_manual(values = colors) +
      xlab("U1") +
      ylab("U2") +
      ggtitle(title) + theme_bw()
  }else{
    colnames(df) = c("U1", "Cluster")
    df$Cluster = as.factor(df$Cluster)
    plots$scat.plot <-
      ggplot(df,
             aes(x = Cluster,
                 y = U1,
                 color = Cluster,
                 fill = Cluster)) +
      geom_violin() +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors) +
      xlab("DBSCAN Clusters") +
      ylab("U1") +
      ggtitle(title) + theme_bw()
  }

  
  # build the KM model
  d <- cbind(time = as.numeric(data_list$time[ind]),
             status = as.numeric(data_list$status[ind]),
             cluster = as.numeric(df$Cluster))
  d <- as.data.frame(d)
  sfit = survfit(Surv(time = d$time,
                      event = d$status) ~cluster,
                 data=d)
  sdiff <- survdiff(Surv(time = d$time,
                         event = d$status) ~cluster,
                    data=d)
  pchsq = pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail = FALSE)
  print(pchsq)
  if(pchsq <= 5e-2){
    #######################
    ### save clustering
    #######################
    barcodes = substr(rownames(embedding[[comboInd]]), start = 9, stop = 12)
    fname = paste0(dbscan_dir, tumortype, "_", comboInd, "_DBSCAN", ".csv")
    cat('\n saving ',fname,'\n')
    write.table(data.frame(barcodes = barcodes,
                           embedding = embedding[[comboInd]],
                           cluster = dbclust$cluster,
                           membership_prob = dbclust$membership_prob,
                           outlier_scores = dbclust$outlier_scores),
                fname)
    
    #######################
    ### save plot of clustering KM
    #######################
    title=paste0("B. Survival Curves")
    plots$surv.plot <-
      autoplot(sfit) +
      scale_colour_manual(values = colors) +
      scale_fill_manual(values = colors) +
      ggtitle(title) +
      theme_bw()
    
    title=text_grob(paste0("TCGA-", tumortype, " Clustering and Survival, UMAP Embedding of:\n ",
                 models$summary$SetEMT[comboInd]," \n ",
                 models$summary$SetINFLAM[comboInd]))
    coxph_cap = text_grob(paste0("A. DBSCAN clustering of ", tumortype, " using gene ontology terms indicative of EMT and INFLAMMATION signatures (min group size = ", mingroupsize, "). \n B. Survival plots corresponding to the clustering on EMT and INFLAMMATION (p = ", format(pchsq, digits=4), ")."))
    fname <- paste0(plotdir, tumortype, "_", comboInd, ".pdf")
    cat('\n saving ',fname,'\n')
    pdf(file = fname, width = 13, height = 7)
    grid.arrange(ggplotGrob(plots$scat.plot),
                 ggplotGrob(plots$surv.plot),
                 ncol = 2,
                 top = title,
                 bottom = coxph_cap)
    dev.off()
    
    #######################
    ### save km latex table
    #######################
    tabl <- data.frame(N = sdiff$n,
                       Observed = sdiff$obs,
                       Expected = sdiff$exp,
                       col5 = (sdiff$obs - sdiff$exp)^2/sdiff$exp,
                       col6 = (sdiff$obs - sdiff$exp)^2/diag(sdiff$var));
    rownames(tabl) = paste0("Cluster ", 1:nrow(sdiff$n));
    tabl = tabl[,-1];
    colnames(tabl) = c("N", "Observed", "Expected", "(O-E)^2/E", "(O-E)^2/V")
    fname = paste0(kmtabdir, tumortype, "_", comboInd, "_km", ".tex")
    captext = paste0("Summary of log-rank test of the Kaplan-Meier estimate for ", tumortype, " on ", models$summary$SetEMT[comboInd]," vs ",models$summary$SetINFLAM[comboInd])
    cat('\n saving ',fname,'\n')
    kmtex = xtable(tabl, label = paste0("table:", tumortype, "_km"), caption = latexify(captext, doublebackslash = F), file = fname)
    print(kmtex, file = fname, compress = FALSE)
  }
  return(pchsq)
}

expandLists <- function(sets,
                        data,
                        cols = ncol(data)){
  crossed_list = data.frame()
  colnames(sets[[2]]) = paste0(colnames(sets[[2]]),'1')
  crossed_list = tidyr::crossing(sets[[1]],sets[[2]])
  names_merged = apply(crossed_list, 1, function(x){
    paste0(x$set_names, ", ", x$set_names1)
  })
  ids_merged = apply(crossed_list, 1, function(x){
    paste0(x$set_ids, ", ", x$set_ids1)
  })
  genes_merged = apply(crossed_list, 1, function(x){
    c(x$genes, x$genes1)
  })
  crossed_list$set_names2 = names_merged
  crossed_list$set_ids2 = ids_merged
  crossed_list$genes2 = genes_merged

  return(crossed_list)
}

embedGenes <- function(data = NULL,
                       data_cols = NULL,
                       genesets = NULL,
                       n_components = NULL){
  # data = the data that is getting embedded, rows are the genes in the various gene sets
  # data_cols = which patients are being embedded
  # genesets = table of gene-set names and gene lists
  embedded = list()
  for(i in 1:length(genesets)){
    ind = match(genesets[[i]], rownames(data))
    raw = data[ind,data_cols]
    set.seed(1)
    embedded[[i]] = umap::umap(t(raw), n_components = n_components)$layout
    print(paste0("finished embedding set ", i, " of ", length(genesets)))
  }
  return(embedded)
}

embedMultiSet <- function(sets = NULL,
                          data = NULL,
                          cols = NULL,
                          tumortype = NULL,
                          outputdir = "../DataRDS/UMAP_embeddings/",
                          n_components = 2){
  if(is.null(cols)){
    cols = 1:ncol(data) # if this is not a fold then use all the data
  }
  expanded = expandLists(sets = sets[[1]],
                              data = data,
                              cols = cols) # create gene lists and name concatenations for every EMT+INFLAM pair
  # get the embedding for EMT
  print(paste0("embedding ", tumortype, " EMT"))
  EMT_embed = embedGenes(data = data,
                         data_cols = cols,
                         geneset = sets[[2]]$EMT$genes,
                         n_components = n_components) # each element is the umap embedding for the i'th EMT set
  # get the embedding for INFLAM
  print(paste0("embedding ", tumortype, " INFLAM"))
  INFLAM_embed = embedGenes(data = data,
                            data_cols = cols,
                            geneset = sets[[3]]$INFLAM$genes,
                            n_components = n_components) # each element is the umap embedding for the i'th INFLAM set
  
  # get the embedding for EMT+INFLAM
  print(paste0("embedding ", tumortype, " EMT+INFLAM"))
  BOTH_embed = embedGenes(data = data,
                          data_cols = cols,
                          geneset = expanded$genes2,
                          n_components = n_components) # each element is the umap embedding for the i'th EMT/INFLAM combo set
  
  embedding = list(EMT = list(embed = EMT_embed,
                              names = sets[[2]]$EMT$set_names,
                              genes = sets[[2]]$EMT$genes),
                   INFLAM = list(embed = INFLAM_embed,
                                 names = sets[[3]]$INFLAM$set_names,
                                 genes = sets[[3]]$INFLAM$genes),
                   BOTH = list(embed = BOTH_embed,
                               names_EMT = expanded$set_names, # each element is the name of the EMT gene set for the i'th EMT+INFLAM combo
                               names_INFLAM = expanded$set_names1, # each element is the name of the INFLAM gene set for the i'th EMT+INFLAM combo
                               genes = expanded$genes2)) 
  saveRDS(embedding, paste0(outputdir, tumortype, "_UMAP_embeddings.RDS"))
  return(embedding)
}