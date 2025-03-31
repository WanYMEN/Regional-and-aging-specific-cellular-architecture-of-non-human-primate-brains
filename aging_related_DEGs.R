## The following pipeline was used for identification of aging-related DEGs ##

suppressMessages({
    library(Seurat)
    library(ggplot2)
})

setwd("/home_path/seurat/harmony")

#### get potential aging-related DEGs across each region and each cell type using FindMarkers ####
## Seurat object after cell-type annotation, which can also be done using the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819; the file included the filtered nucleus versus gene UMI count and data matrices, metadata and UMAP coordinates for submission with a small size ## 
seudt <- readRDS(file = "celltype.harmony.rds")

DefaultAssay(seudt) <- "RNA"

seudt <- SetIdent(seudt, value = "Celltype")

color_ct <- c("SPN" = "#06ACBF", "CGC" = "#F1D5F9", "CB_InN" = "#A643C4", "Bergmann" = "#680886", "TH_ExN" = "#9CF713", "TH_InN" = "#5D9010", "AMY_ExN" = "#147CCB", "CA1-3" = "#D3FF44", "RGC" = "#E105E1", "ExN" = "#FF7F00", "InN" = "#C7960D", "Astro" = "#606FF3", "Micro" = "#5FC45C", "OPC" = "#9E56FD", "cOPC" = "#E2CFFC", "Oligo" = "#F76FB8", "Epen" = "#9B3232", "Endo" = "#A40A0F", "Fib" = "#F65983", "SMC" = "#9A2B48")

cell_types <- names(color_ct)


res_rg_age <- c()
for(i in 1:length(cell_types)){
    seudt_tmp <- subset(seudt, idents = cell_types[i])
    tb_tmp <- data.frame(table(seudt_tmp@meta.data$Age, seudt_tmp@meta.data$Region))
    colnames(tb_tmp) <- c("Age", "Region", "Freq")
    region_tmp <- unique(as.character(tb_tmp[tb_tmp$Freq >=3, ]$Region))

    for(j in 1:length(region_tmp)){
        seudt_tmp2 <- subset(seudt_tmp, subset = Region==region_tmp[j])
        marker_tmp <- c()
        tryCatch({marker_tmp <- FindMarkers(object = seudt_tmp2, ident.1 = "Old", ident.2 = "Young", group.by = "Age", test.use = "MAST", latent.vars = c("nCount_RNA", "CC.Difference", "Sex"), min.pct = 0.25, verbose = FALSE)}, error = function(e) {0})
        if(any(marker_tmp >0)) {
           marker_tmp$gene <- rownames(marker_tmp)
           marker_tmp$celltype <- cell_types[i]
           marker_tmp$region <- region_tmp[j]
           marker_tmp$age1 <- "Old"
           marker_tmp$age2 <- "Young"
           res_rg_age <- rbind(res_rg_age, marker_tmp)
        }
    }
}

write.table(res_rg_age, "allMarkers.age_by_each_celltype_region.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)


old_id <- c("SO1", "SO2")
yound_id <- c("SY1", "SY2")

res_rg_age_id <- c()
for(i in 1:length(cell_types)){
    seudt_tmp <- subset(seudt, idents = cell_types[i])
    for(o in 1:length(old_id)){
        for(y in 1:length(yound_id)){
            seudt_tmp1 <- subset(seudt_tmp, subset = Individual %in% c(old_id[o], yound_id[y]))
            tb_tmp <- data.frame(table(seudt_tmp1@meta.data$Age, seudt_tmp1@meta.data$Region))
            colnames(tb_tmp) <- c("Age", "Region", "Freq")
            region_tmp <- sort(unique(as.character(tb_tmp[tb_tmp$Freq >=3, ]$Region)))

            for(j in 1:length(region_tmp)){
                seudt_tmp2 <- subset(seudt_tmp1, subset = Region==region_tmp[j])
                meta_tmp <- seudt_tmp2@meta.data
                if(nrow(meta_tmp[meta_tmp$Individual==old_id[o],])>=1 & nrow(meta_tmp[meta_tmp$Individual==yound_id[y],])>=1){
                    marker_tmp <- c()
                    tryCatch({marker_tmp <- FindMarkers(object = seudt_tmp2, ident.1 = "Old", ident.2 = "Young", group.by = "Age", test.use = "MAST", latent.vars = c("nCount_RNA", "CC.Difference"), min.pct = 0.25, verbose = FALSE)}, error = function(e) {0})

                    if(any(marker_tmp >0)) {
                       marker_tmp$gene <- rownames(marker_tmp)
                       marker_tmp$celltype <- cell_types[i]
                       marker_tmp$region <- region_tmp[j]
                       marker_tmp$age1 <- old_id[o]
                       marker_tmp$age2 <- yound_id[y]
                       res_rg_age_id <- rbind(res_rg_age_id, marker_tmp)
                    }
                }
            }
        }
    }
}

write.table(res_rg_age_id, "allMarkers.age_by_each_celltype_region_individual.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)



#### filtering to get final candidate aging-related DEGs ####
marker_rg <- read.table("allMarkers.age_by_each_celltype_region.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")

marker_rg.pos <- marker_rg[marker_rg$avg_log2FC >= log2(1.5) & marker_rg$pct.1 >=0.25 & marker_rg$p_val_adj <0.05, ]
marker_rg.neg <- marker_rg[marker_rg$avg_log2FC <= -log2(1.5) & marker_rg$pct.2 >=0.25 & marker_rg$p_val_adj <0.05, ]

marker_rg.pos$DEG <- "pos"
marker_rg.neg$DEG <- "neg"

marker_kept.rg <- unique(rbind(marker_rg.pos, marker_rg.neg))


marker_rg_id <- read.table("allMarkers.age_by_each_celltype_region_individual.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")

marker_rg_id.pos <- marker_rg_id[marker_rg_id$avg_log2FC >= log2(1.5) & marker_rg_id$pct.1 >=0.25 & marker_rg_id$p_val_adj <0.05, ]
marker_rg_id.neg <- marker_rg_id[marker_rg_id$avg_log2FC <= -log2(1.5) & marker_rg_id$pct.2 >=0.25 & marker_rg_id$p_val_adj <0.05, ]

marker_rg_id.pos$DEG <- "pos"
marker_rg_id.neg$DEG <- "neg"

marker_kept.rg_id <- unique(rbind(marker_rg_id.pos, marker_rg_id.neg))

nm_gene.rg_id <- data.frame(ftable(table(marker_kept.rg_id$DEG, marker_kept.rg_id$celltype, marker_kept.rg_id$region, marker_kept.rg_id$gene)))
table(nm_gene.rg_id$Freq)

colnames(nm_gene.rg_id) <- c("DEG", "celltype", "region", "gene", "Freq")
nm_gene.rg_id <- nm_gene.rg_id[nm_gene.rg_id$Freq >=2, ]

final_res <- c()
for(i in 1:nrow(nm_gene.rg_id)){
    if(nm_gene.rg_id[i, ]$Freq==2){
       tmp_rg_id <- marker_kept.rg_id[marker_kept.rg_id$DEG==as.character(nm_gene.rg_id[i, ]$DEG) & marker_kept.rg_id$celltype==as.character(nm_gene.rg_id[i, ]$celltype) & marker_kept.rg_id$region==as.character(nm_gene.rg_id[i, ]$region) & marker_kept.rg_id$gene==as.character(nm_gene.rg_id[i, ]$gene), ]
       
       if(length(unique(tmp_rg_id$age2))==1){
           if(unique(tmp_rg_id$age2)=="SY2"){
              if(unique(tmp_rg_id$region)=="AMY"){
                 tmp_res <- marker_kept.rg[marker_kept.rg$DEG==as.character(nm_gene.rg_id[i, ]$DEG) & marker_kept.rg$celltype==as.character(nm_gene.rg_id[i, ]$celltype) & marker_kept.rg$region==as.character(nm_gene.rg_id[i, ]$region) & marker_kept.rg$gene==as.character(nm_gene.rg_id[i, ]$gene), ]
                 final_res <- rbind(final_res, tmp_res)
              }
           } else {
                tmp_res <- marker_kept.rg[marker_kept.rg$DEG==as.character(nm_gene.rg_id[i, ]$DEG) & marker_kept.rg$celltype==as.character(nm_gene.rg_id[i, ]$celltype) & marker_kept.rg$region==as.character(nm_gene.rg_id[i, ]$region) & marker_kept.rg$gene==as.character(nm_gene.rg_id[i, ]$gene), ]
                final_res <- rbind(final_res, tmp_res)
           }
       } else {
            tmp_res <- marker_kept.rg[marker_kept.rg$DEG==as.character(nm_gene.rg_id[i, ]$DEG) & marker_kept.rg$celltype==as.character(nm_gene.rg_id[i, ]$celltype) & marker_kept.rg$region==as.character(nm_gene.rg_id[i, ]$region) & marker_kept.rg$gene==as.character(nm_gene.rg_id[i, ]$gene), ]
            final_res <- rbind(final_res, tmp_res)
       }
    } else {
       tmp_res <- marker_kept.rg[marker_kept.rg$DEG==as.character(nm_gene.rg_id[i, ]$DEG) & marker_kept.rg$celltype==as.character(nm_gene.rg_id[i, ]$celltype) & marker_kept.rg$region==as.character(nm_gene.rg_id[i, ]$region) & marker_kept.rg$gene==as.character(nm_gene.rg_id[i, ]$gene), ]
       final_res <- rbind(final_res, tmp_res)
   }
}

write.table(final_res, "aging_related_DEGs.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
###### 
