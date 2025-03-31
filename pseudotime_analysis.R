## The following pipeline was used for pseudotime analysis using monocle3, with subtypes mainly from CBC as an example ##

library(igraph") 
library(Seurat)
library(ggplot2)
library(patchwork) 
library(monocle3)
library(dplyr)

setwd("/home_path/monocle3")

data_ct <- data.frame(Celltype = c("CGC", "Purkinje", "Basket/Stellate", "Bergmann"), Label = c("CGC", "Purkinje", "Basket_Stellate", "Bergmann"))
labels <- unique(data_ct$Label)

## Seurat object of nuclei from the cerebellum-specific group using a integrating pipeline similar to 'step3 harmony_integration.R' to identify subtypes##
seudt0 <- readRDS(file = "../seurat/harmony/pos_FindClusters.harmony.CBC.rds")

allmeta <- data.frame(seudt0@meta.data)

color_age <- c("Young" = "#A9A9A9", "Old" = "gray27")

get_earliest_principal_node <- function(cds, ID = "Young"){
    cell_ids <- which(colData(cds)[, "Age"] == ID)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
}


for(i in 1:length(labels)) {
      tmpmeta <- allmeta[(allmeta$Subtype %in% data_ct[data_ct$Label==labels[i],]$Celltype), ]
      if(nrow(tmpmeta) >0){
         seudt <- subset(seudt0, subset = Subtype %in% data_ct[data_ct$Label==labels[i],]$Celltype)
         
         expression_matrix <- seudt@assays$RNA@counts
         metaDat <- data.frame(seudt@meta.data)
         metaDat$Cell <- rownames(metaDat)
         cellmeta <- metaDat[, !grepl("RNA_snn_res", colnames(metaDat))]
         genemeta <- data.frame(gene_short_name = row.names(expression_matrix), row.names = row.names(expression_matrix))
         
         cds <- new_cell_data_set(expression_matrix, cell_metadata = cellmeta, gene_metadata = genemeta)
         if(sum(table(colData(cds)$Age) >10)==2){
            if(table(colData(cds)$Age, colData(cds)$Sex)["Young","female"] >8){
               cds <- preprocess_cds(cds, num_dim = 50)
               cds <- align_cds(cds, alignment_group = "Individual")
               cds <- reduce_dimension(cds)
               cds <- cluster_cells(cds)
               cds <- learn_graph(cds)
               
               pdf(paste("cluster_graph", labels[i], "pdf", sep="."), w=2, h=1.8)
               print(plot_cells(cds, label_principal_points = TRUE, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE))
               dev.off()
               
               cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
               
               saveRDS(cds, file = paste("monocle3", labels[i], "all.rds", sep="."))
            }
         } 
      }
}
#######

