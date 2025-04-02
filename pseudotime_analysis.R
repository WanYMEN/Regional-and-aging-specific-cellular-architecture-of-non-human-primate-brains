## The following pipeline was used for pseudotime analysis using monocle3, with subtypes mainly from CBC as an example ##

library(igraph") 
library(Seurat)
library(ggplot2)
library(patchwork) 
library(monocle3)
library(dplyr)

setwd("/home_path/monocle3")

#### conduct pseudotime analysis using monocle3 for cerebellum-specific subtypes ####
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


#### plot the pseudotemporal trajectory from young towards to aged nuclei of CGC in CBC####
cds <- readRDS(file = "monocle3.CGC.all.rds")

pdf("pseudotime.CGC.time.pdf",w=2.2, h=1.7)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, label_roots = FALSE, graph_label_size=1.5, trajectory_graph_segment_size = 0.3, cell_size = 0.1, rasterize = TRUE) + theme_bw() + theme(panel.background=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.spacing.y = unit(1, 'mm'), legend.text=element_text(size=rel(0.3)), legend.key.size = unit(2, "mm"), legend.title = element_text(size = 5))
dev.off()

pdf("pseudotime.CGC.age.pdf",w=2.2, h=1.7)
plot_cells(cds, label_groups_by_cluster=FALSE, color_cells_by = "Age", label_leaves=FALSE, label_branch_points=FALSE, label_roots = FALSE, trajectory_graph_segment_size = 0.3, cell_size = 0.2) + theme_bw() + theme(panel.background=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.spacing.y = unit(1, 'mm'), legend.text=element_text(size=rel(0.3)), legend.key.size = unit(2, "mm"), legend.title = element_text(size = 5)) + guides(color = guide_legend(override.aes = list(size = 0.9)))
dev.off()


time_dat <- data.frame(Cell = names(cds@principal_graph_aux@listData$UMAP$pseudotime), Pseudotime = as.numeric((cds@principal_graph_aux@listData$UMAP$pseudotime)))
time_dat$Age <- "Young"
time_dat[grep("SO", as.character(time_dat$Cell)), ]$Age <- "Old"
time_dat$Bin <- cut(time_dat$Pseudotime, 20)

num_bin <- as.data.frame.array(table(time_dat$Bin, time_dat$Age))
num_bin$Seq <- 1:20
num_bin$Proportion <- num_bin$Old/(num_bin$Old + num_bin$Young)

pdf("age_ratio_along_pseudotime.CGC.pdf", w=1.2, h=1.2)
ggplot(num_bin, aes(x = Seq, y = Proportion)) + geom_point(size = 0.4) + stat_smooth(method=lm, size = 0.5, se = TRUE) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size = 0.5), axis.text = element_text(size=5, color = "black"), axis.ticks = element_line(colour = "black", size = 0.2), axis.title = element_text(size=6)) + xlab("Pseudotime")
dev.off()

sink("age_ratio_along_pseudotime_lm.CGC.txt")
print(summary(lm(Seq ~ Proportion, data = num_bin)))
sink()

trace('calculateLW', edit = T, where = asNamespace("monocle3")) ## change Matrix::rBind to  rbind, need to run each time
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=1)

pr_deg_ids <- ciliated_cds_pr_test_res %>% filter(q_value < 0.05) %>% arrange(-morans_I)

write.table(pr_deg_ids, "graph_test.CGC.txt", sep="\t", quote=F, col.names = TRUE, row.names = TRUE)

