## The following pipeline was used for identification of cell-type markers ##

suppressMessages({
    library(Seurat)
    library(ggplot2)
})

setwd("/home_path/seurat/harmony")


#### get potential markers using FindAllMarkers ####
options(future.globals.maxSize = 300000 * 1024^2) ## 300GB

## Seurat object after cell-type annotation, which can also be done using the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819; the file included the filtered nucleus versus gene UMI count and data matrices, metadata and UMAP coordinates for submission with a small size ## 
seudt <- readRDS(file = "celltype.harmony.rds")
DefaultAssay(seudt) <- "RNA"
seudt <- SetIdent(seudt, value = "Celltype")


allMarkers <- FindAllMarkers(object = seudt, test.use = "MAST", latent.vars = c("nCount_RNA", "CC.Difference", "Sex"), min.pct = 0.25, only.pos = TRUE, verbose = FALSE)

write.table(allMarkers, "allMarkers.celltype.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)


individuals <- c("SY1", "SY2", "SO1", "SO2")
res_id <- c()
for(i in 1:length(individuals)){
    seudt_tmp <- subset(seudt, subset = Individual==individuals[i])
    marker_tmp <- FindAllMarkers(object = seudt_tmp, test.use = "MAST", latent.vars = c("nCount_RNA", "CC.Difference"), min.pct = 0.25, only.pos = TRUE, verbose = FALSE)
    marker_tmp$individual <- individuals[i]
    res_id <- rbind(res_id, marker_tmp)
}

write.table(res_id, "allMarkers.celltype.byIndividual.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
###### 



#### filtering to get final candidate cell-type markers ####
marker.ct <- read.table("allMarkers.celltype.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
marker_kept.ct <- marker.ct[marker.ct$avg_log2FC >= log2(1.5) & marker.ct$pct.1 >=0.25 & marker.ct$p_val_adj <0.05, ]

marker.ct_id <- read.table("../marker/allMarkers.celltype.byIndividual.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
marker_kept.ct_id <- marker.ct_id[marker.ct_id$avg_log2FC >= log2(1.5) & marker.ct_id$pct.1 >=0.25 & marker.ct_id$p_val_adj <0.05, ]

nm_gene.ct_id <- data.frame(ftable(table(marker_kept.ct_id$cluster, marker_kept.ct_id$gene)))
table(nm_gene.ct_id$Freq)

colnames(nm_gene.ct_id) <- c("cluster", "gene", "Freq")
nm_gene.ct_id <- nm_gene.ct_id[nm_gene.ct_id$Freq >=2, ] ## consider cOPC marker GPR17 only in SO1 and SO2

final_res <- c()
for(i in 1:nrow(nm_gene.ct_id)){
    tmp_res <- marker_kept.ct[marker_kept.ct$cluster==as.character(nm_gene.ct_id[i, ]$cluster) & marker_kept.ct$gene==as.character(nm_gene.ct_id[i, ]$gene), ]
    final_res <- rbind(final_res, tmp_res)
}

write.table(final_res, "celltype_specific_marker.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
###### 
