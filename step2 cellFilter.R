## The following pipeline was used to filter nuclei for each sample ##

suppressMessages({
    library(Matrix)
    library(Seurat)
    library(ggplot2)
    library(DoubletFinder)
    library(getopt)
    })

arg<-matrix(c("input","i","2","character","Path of input",
              "label","l","1","character","Sample label"),
            byrow=T,ncol=5
            )

opt = getopt(arg)

setwd("/home_path/seurat/cellFilter")

## gene annotation information downloaded from BioMart (Ensembl v102) ## 
anno <- read.table("../../data/macaque_v102.bioMart.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")
anno <- anno[, c(1:2,4,9)]
colnames(anno) <- c("EnsID", "GeneName", "Chr", "GeneType")

mtGene <- unique(anno[anno$Chr=="MT",]$EnsID)

## input raw nucleus versus gene UMI count matrix, which can be down from https://doi.org/10.57760/sciencedb.06819 ##
rawdata <- Read10X(data.dir = opt$input, gene.column = 1)

seudt <- CreateSeuratObject(counts = rawdata, project = opt$label, min.cells = 100, min.features = 500)

mt.genes <- intersect(rownames(seudt), mtGene)

seudt[["percent_mt"]] <- PercentageFeatureSet(seudt, features = mt.genes)
seudt[["nCount_nFeature"]] <- seudt[["nCount_RNA"]]/seudt[["nFeature_RNA"]]

seudt <- subset(seudt, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000 & nCount_nFeature >= 1.2 & percent_mt <= 0.1)

seudt <- NormalizeData(seudt, verbose = FALSE)
seudt <- FindVariableFeatures(seudt, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
all.genes <- rownames(seudt)
seudt <- ScaleData(seudt, features = all.genes, vars.to.regress = c("nCount_RNA", "percent_mt"), verbose = FALSE)
seudt <- RunPCA(seudt)
seudt <- RunUMAP(seudt, reduction = "pca", dims = 1:30, verbose = FALSE)
seudt <- RunTSNE(seudt, reduction = "pca", dims = 1:30, perplexity = 30)
seudt <- FindNeighbors(object = seudt, reduction = "pca", dims = 1:30, verbose = FALSE)
seudt <- FindClusters(object = seudt, resolution = 0.5, verbose = FALSE)

## pK Identification (no ground-truth)
sweep_res.seudt <- paramSweep_v3(seudt, PCs = 1:30, sct = FALSE)
sweep_stats.seudt <- summarizeSweep(sweep_res.seudt, GT = FALSE)
bcmvn.seudt <- find.pK(sweep_stats.seudt)

mpK <- as.numeric(as.vector(bcmvn.seudt$pK[which.max(bcmvn.seudt$BCmetric)]))

## Homotypic Doublet Proportion Estimate
annotation.seudt <- seudt@meta.data$seurat_clusters
nExp_poi.seudt <- round(0.1*nrow(seudt@meta.data))  ## Assuming 10% doublet formation rate - tailor for your dataset

## Run DoubletFinder with varying classification stringencies
seudt <- doubletFinder_v3(seudt, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.seudt, reuse.pANN = FALSE, sct = FALSE)

seudt@meta.data$dbnm <- seudt@meta.data[, grep("DF.classifications_", colnames(seudt@meta.data))]

seudt <- subset(seudt, subset = dbnm == "Singlet")

saveRDS(seudt, file = paste("doubletFilter.", opt$label, ".rds", sep = ""))

