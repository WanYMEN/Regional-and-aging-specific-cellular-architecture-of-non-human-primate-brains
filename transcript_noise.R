## The following pipeline was used for transcriptional noise analysis  ##

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(sfsmisc)
    library(MASS)
    library(hopach)
    })

options(future.globals.maxSize = 200000 * 1024^2) ## 200GB

setwd("/home_path/transcript_noise")

## Seurat object after cell-type annotation, which can also be done using the file count_data_meta_UMAP.rds from https://doi.org/10.57760/sciencedb.06819; the file included the filtered nucleus versus gene UMI count and data matrices, metadata and UMAP coordinates for submission with a small size ## 
seudt <- readRDS(file = "../seurat/harmony/celltype.harmony.rds")

color_rg <- c("AMY" = "#ACCAE0", "PU" = "#618D92", "HIP" = "#ADE104", "TH" = "#A9BC8C", "DLPFC" = "#00BFFF", "CG" = "#1A954F", "STG" = "#E1CFEF", "SPL" = "#655CE3", "V4" = "#FF00FF", "CBC" = "#DE7BFC")

regions <- names(color_rg)

color_ct <- c("SPN" = "#06ACBF", "CGC" = "#F1D5F9", "CB_InN" = "#A643C4", "Bergmann" = "#680886", "TH_ExN" = "#9CF713", "TH_InN" = "#5D9010", "AMY_ExN" = "#147CCB", "CA1-3" = "#D3FF44", "V4-ExN" = "#E105E1", "ExN" = "#FF7F00", "InN" = "#C7960D", "Astro" = "#606FF3", "Micro" = "#5FC45C", "OPC" = "#9E56FD", "cOPC" = "#E2CFFC", "Oligo" = "#F76FB8", "Epen" = "#9B3232", "Endo" = "#A40A0F", "Fib" = "#F65983", "SMC" = "#9A2B48")

celltypes <- names(color_ct)

seudt <- SetIdent(seudt, value = "Celltype")
DefaultAssay(seudt) <- "RNA"

all_meta <- seudt@meta.data

### Define function to calculate euclidean distances accounting for cell number and nUMI ####
getEuclideanDistance <- function(celltype, region, lowcv = T){
  print(paste("Working on", celltype, region))
  n_row <- nrow(all_meta[all_meta$Celltype==celltype & all_meta$Region==region, ])
  if(n_row >= 10) {
     tmp <- subset(seudt, subset = Celltype==celltype & Region==region)
     expr <- tmp@assays$RNA@counts
  
     zeros <- which(Matrix::rowSums(expr) == 0)
     expr <- data.matrix(expr[-zeros,])
  
     Down_Sample_Matrix <- function (expr_mat) {
        min_lib_size <- min(colSums(expr_mat))
        down_sample <- function(x) {
           prob <- min_lib_size/sum(x)
           return(unlist(lapply(x, function(y) {
              rbinom(1, y, prob)
           })))
        }
        down_sampled_mat <- apply(expr_mat, 2, down_sample)
        return(down_sampled_mat)
     }
     ds_expr <- Down_Sample_Matrix(expr)
  
     nsample <- min(table(tmp@meta.data$Age))
     if(nsample < 10 | length(table(tmp@meta.data$Age)) < 2){
        print("Not enough cells")
        return(NULL)
     } else {
        old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$Age == "Old")], nsample)
        young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$Age == "Young")], nsample)
        ds_expr_r <- ds_expr[, c(young_r, old_r)]
     
        genes <- c()
        if(lowcv){
           getLowCVgenes <- function(matr){
              means <- Matrix::rowMeans(matr)
              bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
              mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
              asplit <- split(names(means), mean_bin)
              gene_lw <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
                 coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
                 bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
              })))
              return(gene_lw)
           }
           genes <- getLowCVgenes(ds_expr_r)
        } else{
           genes <- rownames(ds_expr_r)
        }
     
        print(length(genes))

        if(length(genes)>0){
           kept_expr <- ds_expr_r[genes, ]  

           calcEuclDist <- function(matr, young, old){
              tmp <- data.matrix(sqrt(matr[, young]))
              mean <- rowMeans(sqrt(matr[, young]))
              d_young <- distancevector(t(tmp), mean , d="euclid")
              names(d_young) <- young
       
              tmp <- data.matrix(sqrt(matr[, old]))
              mean <- rowMeans(sqrt(matr[, old]))
              d_old <- distancevector(t(tmp), mean , d="euclid")
              names(d_old) <- old
       
               return(list(young = d_young, old = d_old))
           }
          ds <- calcEuclDist(matr = kept_expr, old = old_r, young = young_r)
          return(ds)
        } else {
          print("get no LowCV genes")
          return(NULL)
        }
    }
  }
}


res <- lapply(celltypes, function(x) lapply(regions, function(y) getEuclideanDistance(celltype = x, region = y, lowcv = T)))

saveRDS(res, file = "EuclideanDistance.rds")

