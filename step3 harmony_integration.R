## The following pipeline was used to integrate nuclei across the 39 samples ##

suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(plyr)
    library(harmony)
    library(RColorBrewer)
    })

setwd("/home_path/seurat/harmony")

anno <- read.table("../../data/macaque_v102.bioMart.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
anno <- anno[, c(1:2,4,9)]
colnames(anno) <- c("EnsID", "GeneName", "Chr", "GeneType")

mtGene <- unique(anno[anno$Chr=="MT",]$EnsID)

## merged all rds files after DoubletFinder using the merge tool in Seurat
seudt0 <- readRDS(file = "merged.rds")
seudt_count <- GetAssayData(object = seudt0, slot = "counts")
rm(seudt0)

seudt <- CreateSeuratObject(counts = seudt_count, min.cells = 100, min.features = 500)

mt.genes <- intersect(rownames(seudt), mtGene)

seudt[["percent_mt"]] <- PercentageFeatureSet(seudt, features = mt.genes)
seudt[["nCount_nFeature"]] <- seudt[["nCount_RNA"]]/seudt[["nFeature_RNA"]]

seudt <- subset(seudt, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000 & nCount_nFeature >= 1.2 & percent_mt <= 0.1)

region <- unlist(strsplit(rownames(seudt@meta.data), "_"))[3*(1:nrow(seudt@meta.data))-2]
individual <- unlist(strsplit(rownames(seudt@meta.data), "_"))[3*(1:nrow(seudt@meta.data))-1]
sample <- paste(region, individual, sep = "_")

seudt <- AddMetaData(seudt, sample, col.name = "Sample")
seudt <- AddMetaData(seudt, region, col.name = "Region")
seudt <- AddMetaData(seudt, individual, col.name = "Individual")

age <- seudt@meta.data$Individual
age[grepl("SY", age)] <- "Young"
age[grepl("SO", age)] <- "Old"

seudt <- AddMetaData(seudt, age, col.name = "Age")

sex <- seudt@meta.data$Individual
sex[!grepl("SY2", sex)] <- "female"
sex[!grepl("female", sex)] <- "male"

seudt <- AddMetaData(seudt, sex, col.name = "Sex")

## the macaque counterparts of human cell cycle marker genes updated from the list given by Tirosh et al (Tirosh et al. Science 2016, DOI: 10.1126/science.aad0501) ##
s.genes <- read.table("../../data/cellCycle.Tirosh2015.toMc.S.txt", sep="\t", header = FALSE, stringsAsFactors = FALSE, quote="")
g2m.genes <- read.table("../../data/cellCycle.Tirosh2015.toMc.G2M.txt", sep="\t", header = FALSE, stringsAsFactors = FALSE, quote="")

s.genes <- as.character(s.genes$V1)
g2m.genes <- as.character(g2m.genes$V1)

exp_s <-  intersect(s.genes, rownames(seudt))
exp_g2m <-  intersect(g2m.genes, rownames(seudt))

seudt <- NormalizeData(seudt, verbose = FALSE)
seudt <- FindVariableFeatures(seudt, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seudt <- CellCycleScoring(seudt, s.features = exp_s, g2m.features = exp_g2m)
seudt$CC.Difference <- seudt$S.Score - seudt$G2M.Score
all.genes <- rownames(seudt)
seudt <- ScaleData(seudt, features = all.genes, vars.to.regress = c("nCount_RNA", "percent_mt", "CC.Difference", "Sex"), verbose = FALSE)
seudt <- RunPCA(seudt, verbose = FALSE)
seudt <- RunHarmony(seudt, group.by.vars = "Individual")

seudt <- FindNeighbors(seudt, reduction = "harmony", dims = 1:30)
seudt <- FindClusters(seudt, resolution = 1.5, verbose = FALSE)
seudt <- RunUMAP(seudt, reduction = "harmony", dims = 1:30, verbose = FALSE)
seudt <- RunTSNE(seudt, reduction = "harmony", dims = 1:30, perplexity = 50)

saveRDS(seudt, file = "pos_FindClusters.harmony.rds")
