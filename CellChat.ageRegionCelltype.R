## The following pipeline was used for cell-cell communication analysis between cell types across regions for the young and old macaques, respectively ##

suppressMessages({
    library(Seurat)
    library(CellChat)
    library(patchwork)
    })

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 10000 * 1024^2) ## 10GB

setwd("/home_path/CellChat/age")

#### first run CellChat for the young and old groups, respectively ####
## Seurat object after cell-type annotation, which can also be done using the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819; the file included the filtered nucleus versus gene UMI count and data matrices, metadata and UMAP coordinates for submission with a small size ## 
seudt <- readRDS(file = "../../seurat/harmony/celltype.harmony.rds")

DefaultAssay(seudt) <- "RNA"

seudt <- SetIdent(seudt, value = "Celltype")

rg_ct <- paste(seudt@meta.data$Region, seudt@meta.data$Celltype, sep = ".")

seudt <- AddMetaData(seudt, rg_ct, col.name = "RegionCelltype")

seudt_data <- GetAssayData(seudt, slot="data")

all_genes <- rownames(seudt_data)

mc_to_hg <- read.table("/home_path/data/macaque_human.v102.1to1", sep="\t", header = FALSE, stringsAsFactors = FALSE, quote="")
rownames(mc_to_hg) <- mc_to_hg$V1

all_ints <- intersect(all_genes, rownames(mc_to_hg))

seudt_data <- seudt_data[all_ints, ]
rownames(seudt_data) <- mc_to_hg[all_ints, ]$V4

CellChatDB <- CellChatDB.human

CellChatDB.use <- subsetDB(CellChatDB)


seudt_young <- subset(seudt, subset= Age == "Young")
seudt_old <- subset(seudt, subset= Age == "Old")

all_meta_young <- seudt_young@meta.data
cellchat_young <- createCellChat(object = seudt_data[, rownames(all_meta_young)], meta = all_meta_young, group.by = "RegionCelltype")

cellchat_young@DB <- CellChatDB.use

cellchat_young <- subsetData(cellchat_young)

cellchat_young <- identifyOverExpressedGenes(cellchat_young, thresh.pc = 0.25, thresh.fc = log(1.5))
cellchat_young <- identifyOverExpressedInteractions(cellchat_young)

cellchat_young <- computeCommunProb(cellchat_young, population.size = TRUE)
cellchat_young <- filterCommunication(cellchat_young, min.cells = 10)

cellchat_young <- computeCommunProbPathway(cellchat_young)

cellchat_young <- aggregateNet(cellchat_young)

saveRDS(cellchat_young, file = "cellchat_young.rds")


all_meta_old <- seudt_old@meta.data
cellchat_old <- createCellChat(object = seudt_data[, rownames(all_meta_old)], meta = all_meta_old, group.by = "RegionCelltype")

cellchat_old@DB <- CellChatDB.use

cellchat_old <- subsetData(cellchat_old)

cellchat_old <- identifyOverExpressedGenes(cellchat_old, thresh.pc = 0.25, thresh.fc = log(1.5))
cellchat_old <- identifyOverExpressedInteractions(cellchat_old)

cellchat_old <- computeCommunProb(cellchat_old, population.size = TRUE)
cellchat_old <- filterCommunication(cellchat_old, min.cells = 10)

cellchat_old <- computeCommunProbPathway(cellchat_old)

cellchat_old <- aggregateNet(cellchat_old)

saveRDS(cellchat_old, file = "cellchat_old.rds")
########



#### compare cellular communication between the young and old groups ####
color_rg <- c("AMY" = "#ACCAE0", "PU" = "#618D92", "HIP" = "#ADE104", "TH" = "#A9BC8C", "DLPFC" = "#00BFFF", "CG" = "#1A954F", "STG" = "#E1CFEF", "SPL" = "#655CE3", "V4" = "#FF00FF", "CBC" = "#DE7BFC")

regions <- names(color_rg)

color_ct <- c("SPN" = "#06ACBF", "CGC" = "#F1D5F9", "CB_InN" = "#A643C4", "Bergmann" = "#680886", "TH_ExN" = "#9CF713", "TH_InN" = "#5D9010", "AMY_ExN" = "#147CCB", "CA1-3" = "#D3FF44", "RGC" = "#E105E1", "ExN" = "#FF7F00", "InN" = "#C7960D", "Astro" = "#606FF3", "Micro" = "#5FC45C", "OPC" = "#9E56FD", "cOPC" = "#E2CFFC", "Oligo" = "#F76FB8", "Epen" = "#9B3232", "Endo" = "#A40A0F", "Fib" = "#F65983", "SMC" = "#9A2B48")

celltypes <- names(color_ct)

color_age <- c("Young" = "#A9A9A9", "Old" = "gray27")

all_rg_ct <- paste(rep(regions, each=length(celltypes)), celltypes, sep = ".")

all_colors <- distinctColorPalette(k = length(all_rg_ct), altCol = FALSE, runTsne = FALSE)
names(all_colors) <- all_rg_ct

cellchat_young <- readRDS(file = "cellchat_young.rds")
cellchat_old <- readRDS(file = "cellchat_old.rds")

object.list <- list(Young = cellchat_young, Old = cellchat_old)

for (i in 1:length(object.list)) {
    object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
}

group.new <- union(levels(object.list[[1]]@idents), levels(object.list[[2]]@idents))
group.new <- all_rg_ct[all_rg_ct %in% group.new]
object.list[[1]] <- liftCellChat(object.list[[1]], group.new)
object.list[[2]] <- liftCellChat(object.list[[2]], group.new)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pdf("compareInteractions.age.pdf", w=3.5, h=2.5)
compint_count <- compareInteractions(cellchat, show.legend = F, group = c(1,2), size.text = 7)
compint_weight <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight", size.text = 7)
compint_count + compint_weight
dev.off()


options(ggrepel.max.overlaps = 30)

pdf("diff_signalingRole_scatter.age.pdf", w=2.6, h=2.6)
netAnalysis_diff_signalingRole_scatter(cellchat, comparison = c(1, 2), color.use = all_colors[group.new], label.size = 2, font.size = 7, font.size.title = 9)
dev.off()

## Compare the overall information flow of each signaling pathway or ligand-receptor pair
pdf("rankNet_bar.age.pdf", w=8.3, h=2)
rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, font.size = 6, tol = 0.1, do.flip = FALSE, x.angle = 90)
dev.off()
########

