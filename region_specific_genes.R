## The following pipeline was used for identification of region-specific genes ##

suppressMessages({
    library(Seurat)
    library(ggplot2)
})

setwd("/home_path/seurat/harmony")

#### get potential markers of a region for each cell type using FindMarkers ####
## Seurat object after cell-type annotation, which can also be done using the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819; the file included the filtered nucleus versus gene UMI count and data matrices, metadata and UMAP coordinates for submission with a small size ## 
seudt <- readRDS(file = "celltype.harmony.rds")

DefaultAssay(seudt) <- "RNA"

seudt <- SetIdent(seudt, value = "Celltype")

color_ct <- c("SPN" = "#06ACBF", "CGC" = "#F1D5F9", "CB_InN" = "#A643C4", "Bergmann" = "#680886", "TH_ExN" = "#9CF713", "TH_InN" = "#5D9010", "AMY_ExN" = "#147CCB", "CA1-3" = "#D3FF44", "V4-ExN" = "#E105E1", "ExN" = "#FF7F00", "InN" = "#C7960D", "Astro" = "#606FF3", "Micro" = "#5FC45C", "OPC" = "#9E56FD", "cOPC" = "#E2CFFC", "Oligo" = "#F76FB8", "Epen" = "#9B3232", "Endo" = "#A40A0F", "Fib" = "#F65983", "SMC" = "#9A2B48")

celltypes <- names(color_ct)


res_rg <- c()
for(i in 1:length(celltypes)){
    seudt_tmp <- subset(seudt, idents = celltypes[i])
    tb_tmp <- data.frame(table(seudt_tmp@meta.data$Region))
    colnames(tb_tmp) <- c("Region", "Freq")
    region_tmp <- unique(as.character(tb_tmp[tb_tmp$Freq >=3, ]$Region))

    for(j in 1:length(region_tmp)){
        marker_tmp <- FindMarkers(object = seudt_tmp, ident.1 = region_tmp[j], group.by = "Region", test.use = "MAST", latent.vars = c("nCount_RNA", "CC.Difference", "Sex"), min.pct = 0.25, only.pos = TRUE, verbose = FALSE)
        marker_tmp$gene <- rownames(marker_tmp)
        marker_tmp$celltype <- celltypes[i]
        marker_tmp$region <- region_tmp[j]
        res_rg <- rbind(res_rg, marker_tmp)
    }
}

write.table(res_rg, "allMarkers.region_by_each_celltype.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)


individuals <- c("SY1", "SY2", "SO1", "SO2")

res_rg_id <- c()
for(i in 1:length(celltypes)){
   for(k in 1:length(individuals)){
       seudt_tmp <- subset(seudt, idents = celltypes[i])
       seudt_tmp <- subset(seudt_tmp, subset = Individual==individuals[k])
       tb_tmp <- data.frame(table(seudt_tmp@meta.data$Region))
       colnames(tb_tmp) <- c("Region", "Freq")
       region_tmp <- unique(as.character(tb_tmp[tb_tmp$Freq >=3, ]$Region))

       for(j in 1:length(region_tmp)){
           ot_tmp <- seudt_tmp@meta.data[seudt_tmp@meta.data$Region != region_tmp[j],]

           if(nrow(ot_tmp) >=3){
              marker_tmp <- c()
              tryCatch({marker_tmp <- FindMarkers(object = seudt_tmp, ident.1 = region_tmp[j], group.by = "Region", min.cells.group = 1, test.use = "MAST", latent.vars = c("nCount_RNA", "CC.Difference"), min.pct = 0.25, only.pos = TRUE, verbose = FALSE)}, error = function(e) {0})
              if(any(marker_tmp >0)) {
                 marker_tmp$gene <- rownames(marker_tmp)
                 marker_tmp$celltype <- celltypes[i]
                 marker_tmp$region <- region_tmp[j]
                 marker_tmp$individual <- individuals[k]
                 res_rg_id <- rbind(res_rg_id, marker_tmp)
              }
           }
       }
  }
}

write.table(res_rg_id, "allMarkers.region_by_each_celltype_individual.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
###### 



#### filtering to get final candidate region-specific genes ####
regions <- c("AMY", "PU", "HIP", "TH", "DLPFC", "CG", "STG", "SPL", "V4", "CBC")

marker.rg_ct <- read.table("allMarkers.region_by_each_celltype.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
marker_kept.rg_ct <- marker.rg_ct[marker.rg_ct$avg_log2FC >= log2(1.5) & marker.rg_ct$pct.1 >=0.25 & marker.rg_ct$p_val_adj <0.05, ]

kept_rg_ct <- c()
for(i in 1:length(celltypes)){
    rg_eachCT <- marker_kept.rg_ct[marker_kept.rg_ct$celltype==celltypes[i], ]
    if(nrow(rg_eachCT) >0){
        dup_ct <- rg_eachCT[duplicated(rg_eachCT$gene), ]$gene
        rg_eachCT <- rg_eachCT[!(rg_eachCT$gene %in% dup_ct), ]
        kept_rg_ct <- rbind(kept_rg_ct, rg_eachCT)
    }
}


marker.rg_ct_id <- read.table("allMarkers.region_by_each_celltype_individual.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
marker_kept.rg_ct_id <- marker.rg_ct_id[marker.rg_ct_id$avg_log2FC >= log2(1.5) & marker.rg_ct_id$pct.1 >=0.25 & marker.rg_ct_id$p_val_adj <0.05, ]
nm_gene.rg_ct_id <- data.frame(ftable(table(marker_kept.rg_ct_id$celltype, marker_kept.rg_ct_id$region, marker_kept.rg_ct_id$gene)))
table(nm_gene.rg_ct_id$Freq)

colnames(nm_gene.rg_ct_id) <- c("celltype", "region", "gene", "Freq")
nm_gene.rg_ct_id <- nm_gene.rg_ct_id[nm_gene.rg_ct_id$Freq >=2, ]

sp_rg <- c()
for(i in 1:nrow(nm_gene.rg_ct_id)){
    tmp_res <- kept_rg_ct[kept_rg_ct$celltype==as.character(nm_gene.rg_ct_id[i, ]$celltype) & kept_rg_ct$region==as.character(nm_gene.rg_ct_id[i, ]$region) & kept_rg_ct$gene==as.character(nm_gene.rg_ct_id[i, ]$gene), ]
    
    sp_rg <- rbind(sp_rg, tmp_res)
}


sp_ct <- read.table("celltype_specific_marker.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
sp_ct$celltype <- sp_ct$cluster
sp_ct <- sp_ct[, -6]

AMY <- c("AMY_ExN")
CBC <- c("CGC", "CB_InN", "Bergmann")
HIP <- c("CA1-3")
V4 <- c("V4-ExN")
PU <- c("SPN")
TH <- c("TH_ExN", "TH_InN")

all_list <- list(AMY, CBC, HIP, V4, PU, TH)
names(all_list) <- c("AMY", "CBC", "HIP", "V4", "PU", "TH")

res_ct <- c()
for(i in 1:length(all_list)){
    tmp_ct_tg <- sp_ct[sp_ct$celltype %in% all_list[[i]], ]
    tmp_ct_ot <- sp_ct[!(sp_ct$celltype %in% all_list[[i]]), ]
    tmp_ct_tg <- tmp_ct_tg[!(tmp_ct_tg$gene %in% unique(tmp_ct_ot$gene)), ]
    tmp_ct_tg$region <- names(all_list)[i]
    res_ct <- rbind(res_ct, tmp_ct_tg)
}

res_ct$category <- paste(res_ct$celltype, res_ct$region, res_ct$gene, sep = "_")
sp_rg$category <- paste(sp_rg$celltype, sp_rg$region, sp_rg$gene, sep = "_")

res_ct <- res_ct[!(res_ct$category %in% unique(sp_rg$category)), ]
all_res <- rbind(sp_rg, res_ct)

write.table(all_res, "region_specific_marker.txt", sep="\t", quote=F, col.names = TRUE, row.names = FALSE)
###### 
