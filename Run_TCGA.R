remove(list = ls())
source("TCGAbiolinks_all_function.R")
library(gridExtra)
par(mfrow=c(2,2))
n_clust <- 2
df <- n_clust - 1
data <- gdcQ()

## Cluster and Plot the EMT/Inflam axes
kw1 <- list(EMT = c("EPITHELIAL_TO_MESENCHYMAL",
                    "EMT"),
            INFLAM = c("INFLAMMATORY",
                       "INFLAMMATION"))
sets1 <- kwsearch(keywords = kw1,
                  data_genes = data$data_genes)
pcs1 <- findPC(data,sets1)
clust1 <- clusterPC(data,pcs1,n_clust)
surv1 <- survMax(clust1,
                data$days_to_death,
                data$vital_status,
                deg_freedom = df)
p_both <- 
  survivalPlotFig(data$vital_status,
                data$days_to_death,
                clust1,
                pcs1,
                surv1$min_index,
                axes = c(1,2),
                violin = FALSE,
                cols = c("#a6cee3",
                         "#1f78b4"))



# Cluster and Plot EMT only
kw2 <- list(EMT = c("EPITHELIAL_TO_MESENCHYMAL","EMT"))
sets2 <- kwsearch(keywords = kw2,
                 data_genes = data$data_genes)
pcs2 <- findPC(data,sets2)
clust2 <- clusterPC(data,pcs2,n_clust)
surv2 <- survMax(clust2,
                data$days_to_death,
                data$vital_status,
                deg_freedom = df)
p_emt <- 
  survivalPlotFig(data$vital_status,
                data$days_to_death,
                clust2,
                pcs2,
                surv2$min_index,
                axes = c(1),
                violin = TRUE,
                cols = c("#b2df8a",
                         "#33a02c"))

grid.arrange(p_both$scat.plot, 
             p_both$surv.plot$plot, 
             p_emt$scat.plot, 
             p_emt$surv.plot$plot, 
             ncol = 2,
             top = "Clustering and Survival")

