
library(ggplot2)

setwd("/home_path/seurat/harmony")

## gene length information downloaded from Ensembl (v102) ## 
gene_length <- read.table(file = "../../data/Macaca_mulatta.Mmul_10.v102.gene_length.txt", header = TRUE, row.names = 1)

## spatial distance values between any two of the ten brain regions ## 
dis_rg <- read.table("../../data/region_distance.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="")

## gene versus summed UMI count matrix per region and individual  ## 
data_pseudo <- read.table("../../data/count_sum_region_individual.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE, quote="", row.names = 1)

gene_length <- gene_length[rownames(data_pseudo), ]
total_counts <- colSums(x = data_pseudo)
gene_length_kb <- as.vector(x = gene_length / 1000)
pseudo_fpkm <- sweep(data_pseudo, 2, total_counts, FUN = "/") * 1e6 / gene_length_kb

kept_bulk <- pseudo_fpkm[apply(pseudo_fpkm, 1, function(x) sum(x>=1) >= 2), ]
exp_bulk <- log2(kept_bulk + 1)

individuals <- c("SY1", "SY2", "SO1", "SO2")
regions <- c("AMY", "PU", "HIP", "TH", "DLPFC", "CG", "STG", "SPL", "V4", "CBC")

data_dis <- c()
for(i in 1:length(individuals)){
    for(j in 1:(length(regions)-1)){
        for(k in (j+1):length(regions)){
            samp1 <- paste(regions[j], individuals[i], sep = "_")
            samp2 <- paste(regions[k], individuals[i], sep = "_")
            
            if(samp1!="AMY_SY1" & samp2!="AMY_SY1"){
               tmp_exp_dis <- dist(rbind(as.numeric(exp_bulk[, samp1]), as.numeric(exp_bulk[, samp2])), method = "euclidean")
               tmp_res <- c(individuals[i], regions[j], regions[k], tmp_exp_dis)
               data_dis <- rbind(data_dis, tmp_res)
            }
        }
    }
}

data_dis <- data.frame(data_dis)
colnames(data_dis) <- c("Individual", "Region1", "Region2", "Dis_exp")
data_dis$Dis_exp <- as.numeric(data_dis$Dis_exp)

data_dis$Dis_rg <- NA
for(i in 1:nrow(data_dis)){
    data_dis[i,]$Dis_rg <- dis_rg[as.character(data_dis[i, ]$Region1), as.character(data_dis[i, ]$Region2)]
}

data_dis$Age <- "Young"
data_dis[grep("SO", data_dis$Individual), ]$Age <- "Old"

color_age <- c("Young" = "#A9A9A9", "Old" = "gray27")

pdf("dis_exp.pdf", w=1.6, h=1.8)
ggplot(data_dis, aes(x = Dis_rg, y = Dis_exp)) + geom_point(size = 0.6) + stat_smooth(method=lm, size = 0.6, se = TRUE) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size = 0.6), axis.text = element_text(size=6, color = "black"), axis.ticks = element_line(colour = "black", size = 0.4), axis.title = element_text(size=7)) + xlab("Region distance") + ylab("Expression distance")
dev.off()

pdf("dis_exp.age.pdf", w=3, h=1.9)
ggplot(data_dis, aes(x = Dis_rg, y = Dis_exp)) + geom_point(size = 0.6, aes(color = Age)) + stat_smooth(method=lm, size = 0.6, se = FALSE, aes(color = Age, group = Age)) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size = 0.6), axis.text = element_text(size=6, color = "black"), axis.ticks = element_line(colour = "black", size = 0.6), axis.title.x = element_text(size=7), axis.title.y = element_text(size=9)) + xlab("Region distance") + ylab("Expression distance") + scale_colour_manual(values = color_age)
dev.off()

sink("dis_exp.lm.txt")
print("All")
print(summary(lm(Dis_exp ~ Dis_rg, data = data_dis)))
print("")
print("Young")
print(summary(lm(Dis_exp ~ Dis_rg, data = data_dis[data_dis$Age=="Young", ])))
print("")
print("Old")
print(summary(lm(Dis_exp ~ Dis_rg, data = data_dis[data_dis$Age=="Old", ])))
sink()

