## The following pipeline was used for the UMAP plot of cell types in the Fig. 1B ##

suppressMessages({
    library(Seurat)
    library(ggplot2)
})

setwd("/home_path/seurat/harmony")

## Seurat object after cell-type annotation, which can also be done using the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819; the file included the filtered nucleus versus gene UMI count and data matrices, metadata and UMAP coordinates for submission with a small size ## 
seudt <- readRDS(file = "celltype.harmony.rds")

color_ct_subtype <- c("SPN" = "#06ACBF", "CGC" = "#F1D5F9", "Purkinje" = "#D157F6", "Basket/Stellate" = "#9D56B3", "Bergmann" = "#680886", "TH_ExN" = "#9CF713", "TH_InN" = "#5D9010", "AMY_ExN" = "#147CCB", "CA1-3" = "#D3FF44", "V4-ExN" = "#E105E1", "L2/3" = "#FCF2E8", "L2/3/4" = "#FFE8D1", "L4" = "#FCC58F", "UL IT" = "#FFA64D", "L5" = "#EF7B07", "L5/6" = "#CB843D", "L6" = "#B5610D", "L5/6 NP" = "#98612A", "L6b" = "#864505", "DL CT" = "#653F1A", "InN_SST" = "#FAFA36", "InN_PVALB" = "#C6AC08", "InN_ADARB2" = "#C0B085", "InN_VIP" = "#44433A", "InN_LAMP5" = "#FAC533", "Astro" = "#606FF3", "Micro" = "#5FC45C", "OPC" = "#9E56FD", "cOPC" = "#E2CFFC", "Oligo" = "#F76FB8", "Epen" = "#9B3232", "Endo" = "#A40A0F", "Fib" = "#F65983", "SMC" = "#9A2B48")

DefaultAssay(seudt) <- "RNA"

seudt@meta.data$Subtype <- factor(seudt@meta.data$Subtype, levels = names(color_ct_subtype), order = TRUE)

seudt <- SetIdent(seudt, value = "Subtype")

pdf("UMAP.subtype.pdf", w=5.4, h=3.4)
DimPlot(seudt, reduction = "umap", label = TRUE, label.size = 1.8, cols = color_ct_subtype, order = TRUE, repel = TRUE) + theme(legend.text=element_text(size=rel(0.6)), legend.key.size = unit(1, "mm"), legend.spacing.y = unit(0.1, 'mm')) + guides(color = guide_legend(override.aes = list(size=1.5), ncol=2))
dev.off()
