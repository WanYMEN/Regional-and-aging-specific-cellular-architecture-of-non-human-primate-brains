## The following pipeline was used for cell-cell communication analysis between cell types ##

suppressMessages({
    library(Seurat)
    library(CellChat)
    library(patchwork)
    })

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 300000 * 1024^2) ## 300GB

setwd("/home_path/CellChat/celltype")

## Seurat object after cell-type annotation, which can also be done using the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819; the file included the filtered nucleus versus gene UMI count and data matrices, metadata and UMAP coordinates for submission with a small size ## 
seudt <- readRDS(file = "../../seurat/harmony/celltype.harmony.rds")

DefaultAssay(seudt) <- "RNA"

seudt <- SetIdent(seudt, value = "Celltype")

seudt_data <- GetAssayData(seudt, slot="data")

all_genes <- rownames(seudt_data)

## One-to-one orthologous genes between human and rhesus macaque were identified using BioMart (Ensembl v102) ##
mc_to_hg <- read.table("/home_path/data/macaque_human.v102.1to1", sep="\t", header = FALSE, stringsAsFactors = FALSE, quote="")
rownames(mc_to_hg) <- mc_to_hg$V1

all_ints <- intersect(all_genes, rownames(mc_to_hg))

seudt_data <- seudt_data[all_ints, ]
rownames(seudt_data) <- mc_to_hg[all_ints, ]$V4

all_meta <- seudt@meta.data

cellchat <- createCellChat(object = seudt_data, meta = all_meta, group.by = "Celltype")

groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database (default CellChatDB) in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost. This step is necessary even if using the whole database
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 6) # do parallel, cannot set too large, or need lots of memory

cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = 0.25, thresh.fc = log(1.5))
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, population.size = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "cellchat_celltype.rds")

