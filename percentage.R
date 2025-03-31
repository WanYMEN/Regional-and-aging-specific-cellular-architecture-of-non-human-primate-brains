## The following pipeline was used to calculate the percentages of nuclei per cell type across regions, age groups or individuals, respectively ##

library(ggplot2)
library(reshape)
library(dplyr)

setwd("/home_path/seurat/harmony")

## metadata information, which can also be obtained from the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819 ## 
metaDt <- read.table("metadata.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")

median(metaDt$nCount_RNA) 
median(metaDt$nFeature_RNA)
mean(metaDt$nFeature_RNA) 

individuals <- c("SY1", "SY2", "SO1", "SO2")
regions <- c("AMY", "PU", "HIP", "TH", "DLPFC", "CG", "STG", "SPL", "V4", "CBC")

color_ct <- c("SPN" = "#06ACBF", "CGC" = "#F1D5F9", "CB_InN" = "#A643C4", "Bergmann" = "#680886", "TH_ExN" = "#9CF713", "TH_InN" = "#5D9010", "AMY_ExN" = "#147CCB", "CA1-3" = "#D3FF44", "RGC" = "#E105E1", "ExN" = "#FF7F00", "InN" = "#C7960D", "Astro" = "#606FF3", "Micro" = "#5FC45C", "OPC" = "#9E56FD", "cOPC" = "#E2CFFC", "Oligo" = "#F76FB8", "Epen" = "#9B3232", "Endo" = "#A40A0F", "Fib" = "#F65983", "SMC" = "#9A2B48")

celltypes <- names(color_ct)

color_age <- c("Young" = "#A9A9A9", "Old" = "gray27")
color_rg <- c("AMY" = "#ACCAE0", "PU" = "#618D92", "HIP" = "#ADE104", "TH" = "#A9BC8C", "DLPFC" = "#00BFFF", "CG" = "#1A954F", "STG" = "#E1CFEF", "SPL" = "#655CE3", "V4" = "#FF00FF", "CBC" = "#DE7BFC")
color_id <- c("SY1" = "#BDB76B", "SY2" = "#808000", "SO1" = "#FF3399", "SO2" = "#8B008B")


## get nuclei number of each cell type ## 
num_celltype <- as.data.frame(table(metaDt$Celltype))
colnames(num_celltype) <- c("Celltype", "Count")
num_celltype$Celltype <- factor(num_celltype$Celltype, levels = rev(celltypes), order = TRUE)


## the percentages of nuclei per cell type across individuals ## 
num_individual <- data.frame(table(metaDt[, c("Celltype", "Individual", "Age")]))
colnames(num_individual) <- c("Celltype", "Individual", "Age", "Count")
num_individual <- num_individual[num_individual$Count >0,]
mt_individual <- num_individual %>% group_by(Individual) %>% mutate(Percentage0 = prop.table(Count)) %>% data.frame
mt_individual_age <- mt_individual %>% group_by(Age) %>% mutate(Percentage1 = prop.table(Percentage0)) %>% data.frame
mt_individual_ct <- mt_individual_age %>% group_by(Celltype) %>% mutate(Percentage = prop.table(Percentage1)) %>% data.frame
mt_individual_ct$Celltype <- factor(mt_individual_ct$Celltype, levels = rev(celltypes), order = TRUE)
mt_individual_ct$Individual <- factor(mt_individual_ct$Individual, levels = individuals, order = TRUE)


## the percentages of nuclei per cell type across regions ## 
num_region <- data.frame(table(metaDt[, c("Celltype", "Region")]))
colnames(num_region) <- c("Celltype", "Region", "Count")
num_region <- num_region[num_region$Count >0,]
mt_region <- num_region %>% group_by(Region) %>% mutate(Percentage0 = prop.table(Count)) %>% data.frame
mt_region_ct <- mt_region %>% group_by(Celltype) %>% mutate(Percentage = prop.table(Percentage0)) %>% data.frame
mt_region_ct$Celltype <- factor(mt_region_ct$Celltype, levels = rev(celltypes), order = TRUE)
mt_region_ct$Region <- factor(mt_region_ct$Region, levels = regions, order = TRUE)


## the percentages of nuclei per cell type across two age groups ## 
num_age <- data.frame(table(metaDt[, c("Celltype", "Age")]))
colnames(num_age) <- c("Celltype", "Age", "Count")
num_age <- num_age[num_age$Count >0,]
mt_age <- num_age %>% group_by(Age) %>% mutate(Percentage0 = prop.table(Count)) %>% data.frame
mt_age_ct <- mt_age %>% group_by(Celltype) %>% mutate(Percentage = prop.table(Percentage0)) %>% data.frame
mt_age_ct$Celltype <- factor(mt_age_ct$Celltype, levels = rev(celltypes), order = TRUE)
mt_age_ct$Age <- factor(mt_age_ct$Age, levels = names(color_age), order = TRUE)

