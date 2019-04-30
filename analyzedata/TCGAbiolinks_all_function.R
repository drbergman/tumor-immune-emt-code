library(SummarizedExperiment)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(ggplot2)
library(msigdbr)
library(tidyr)
library(gridExtra)


gdcQ <- function(proj = "TCGA-PAAD"){
  ## get the expression and survival data
  query <- GDCquery(project = proj,
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  GDCdownload(query, method = "api", files.per.chunk = 10)
  data <- GDCprepare(query)
  
  data_matrix <- assay(data, 'normalized_count')
  data_genes <- rownames(data_matrix)
  vital_status <- data$vital_status
  days_to_death <- data$days_to_death
  # assume the data is censored (all alive patients follow up by the time the last one dies)
  for(i in 1:length(days_to_death)){
    if(is.na(days_to_death[i])){
      days_to_death[i] <- max(sort(days_to_death))
    }
  }
  
  # binarize the vital status
  for(i in 1:length(vital_status)){
    if(vital_status[i] == 'alive'){
      vital_status[i] <- 0
    } else {
      vital_status[i] <- 1
    }
  }
  
  return(list(data = data_matrix,
              data_genes = data_genes,
              vital_status = vital_status,
              days_to_death = days_to_death))
}

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

## project the data on per-gene-list components (PC1 from each respective dim reduction)
## creates a list of data frames, one data frame per kw set
## col 1 is the gene set name, col2 is PC1 for pca on the data across genes in the gene set

findPC <- function(data_list,
                       set_list){
  #browser()
  ## make a nested list
  ##  level 1: keyword groups
  ##  level 2: PC1 of each gene set for that keyword
  pc_list <- list()
  for(i in 1:length(set_list)){
    gs_table <- data.frame(matrix(0, nrow(set_list[[i]]),0))
    pc_data <- list()
    for(j in 1:nrow(set_list[[i]])){
      rows <- unlist(set_list[[i]]$genes[j])
      data_gene_set <- data_list$data[rows,]
      pc_data[[j]] <- 
        prcomp(data_gene_set)$rotation[,1]
    }
    gs_table$set_names <- set_list[[i]]$set_names
    gs_table$set_ids <- set_list[[i]]$set_ids
    gs_table$pc_data <- pc_data
    pc_list[[as.name(names(set_list)[i])]] <- gs_table
  }
  return(pc_list)
}

## cluster the patients based on n-dimensional data 
## look at all possible combinations of gene sets from the 
## keyword groups
## calculate the number of 

clusterPC <- function(data_list,
                      pc_list,
                      n_clusters){
  #browser()
  n_groups <- length(pc_list)
  product_table <- data.frame()
  for(i in 1:n_groups){
    product_table <- tidyr::crossing(product_table, pc_list[[i]])
  }
  
  cluster_table <- data.frame(matrix(0,nrow(product_table),0))
  names <- list()
  ids <- list()
  clusters <- list()
  name_idx <- seq(1, ncol(product_table), 3)
  id_idx <- seq(2, ncol(product_table), 3)
  pc_idx <- seq(3, ncol(product_table), 3)
  for(i in 1:nrow(product_table)){
    names[[i]] <- unlist(product_table[i,name_idx])
    names(names[[i]]) <- names(pc_list)
    ids[[i]] <- unlist(product_table[i,id_idx])
    names(ids[[i]]) <- names(pc_list)

    pc <- product_table[i,pc_idx]
    clusters[[i]] <- kmeans(cbind(sapply(pc, unlist)),n_clusters)$cluster
  }
  cluster_table$names <- names
  cluster_table$ids <- ids
  cluster_table$clusters <- clusters
  return(cluster_table)
}

survMax <- function(cluster_table,
                    days_to_death,
                    vital_status,
                    deg_freedom){
  #browser()
  surv_list <- list()
  for(i in 1:nrow(cluster_table)){
    d <- cbind(days_to_death = as.numeric(days_to_death),
               vital_status = as.numeric(vital_status),
               
               
               
               #### need fix the col index ####
               
               
               cluster = cluster_table$clusters[[i]])
    sfit <- survfit(Surv(time = days_to_death,
                         event = vital_status) ~cluster,
                    data=as.data.frame(d))
    sdiff <- survdiff(Surv(time = days_to_death,
                           event = vital_status) ~cluster,
                      data=as.data.frame(d))
    pval <- pchisq(sdiff$chisq,
                   df = deg_freedom,
                   lower = FALSE)
    
    surv_list[[i]] <- list(sfit = sfit,
                           sdiff = sdiff,
                           plogrank = pval)
  }
  
  # get the min p-value pair
  min <- 1e20
  optimum <- 0
  for(i in 1:length(surv_list)){
    if(surv_list[[i]]$plogrank < min){
      min <- surv_list[[i]]$plogrank
      optimum <- i
    }
  }
  
  return(list(surv_list = surv_list,
              min_index = optimum,
              min_p = min))
}

# makes a survival plot and cluster plot for a single kaplan meier model
survivalPlotFig <- function(vital_status,
                            days_to_death,
                            cluster_table,
                            pc_list,
                            index,
                            axes = c(1,2),
                            violin = FALSE,
                            cols){
  # plot the clusters and survival
  palette(value=cols)
  #browser()
  sets <- cluster_table$names[[index]][axes] # only two dimensions may be plotted
  set_ids <- cluster_table$ids[[index]][axes]
  n_cases <- length(cluster_table$clusters[[1]])
  pc_table <- data.frame(matrix(0, n_cases, 0))
  for(i in 1:length(sets)){
    group <- pc_list[[i]]$set_names
    set_idx <- which(group == sets[i])
    pc <- pc_list[[i]]$pc_data[set_idx]
    pc_table <- cbind(pc_table, mynewvar = pc)
    colnames(pc_table)[ncol(pc_table)] <- set_ids[i]
  }
  
  pc_table <- cbind(pc_table, cluster = cluster_table$clusters[[index]])
  pc_table$cluster <- as.factor(pc_table$cluster)
  if(violin){
    scat.plot <-
      ggplot(pc_table,
           aes_(x = ~`cluster`,
                y = as.name(colnames(pc_table)[1]),
                fill = ~`cluster`)) +
      geom_violin(trim=FALSE) +
      xlab("cluster") +
      ylab(paste0(names(set_ids)[1], ": ", set_ids[1])) +
      geom_jitter(shape=16, position=position_jitter(0.2)) +
      scale_fill_manual(values=cols) +
      ggtitle(paste0('Clustering on ',
                     names(set_ids)[1])) +
      theme_bw()
  }else{
    scat.plot <-
      ggplot(pc_table,
             aes_(x = as.name(colnames(pc_table)[1]),
                  y = as.name(colnames(pc_table)[2]),
                  color = ~`cluster`)) +
        geom_point() +
        xlab(paste0(names(set_ids)[1], ": ", set_ids[1])) +
        ylab(paste0(names(set_ids)[2], ": ", set_ids[2])) +
        scale_color_manual(values = cols) +
        ggtitle(paste0('Clustering on ',
                       names(set_ids)[1],
                       ' and ',
                       names(set_ids)[2])) +
        theme_bw()
  }
  d <- cbind(days_to_death = as.numeric(days_to_death),
             vital_status = as.numeric(vital_status),
             cluster = as.numeric(pc_table$cluster))
  survdata <- as.data.frame(d)
  sfit <- survfit(Surv(time = days_to_death,
                       event = vital_status) ~cluster,
                  data = survdata)
  if(violin){
    title=paste0("Survival Curves for cluster on \n",
                 names(set_ids)[1])
  }else{
    title=paste0("Survival Curves for cluster on \n", 
                 names(set_ids)[1],
                 " and ",
                 names(set_ids)[2])
  }
  surv.plot <-
    ggsurvplot(sfit,
             title=title,
             palette=cols,
             legend.labs=c("cluster 1", "cluster 2"),
             data = survdata,
             xlab = "Days")
  return(list(scat.plot = scat.plot,
              surv.plot = surv.plot))
  
}