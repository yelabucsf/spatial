# load packages
library("reticulate")
use_python("/Applications/anaconda3/bin/python3", required = T)
sc <- import("scanpy")
library(tidyverse)
library(NMF)
install.extras('NMF')

# load liver data (decontaminated and preprocessed)
setwd("~/Documents/projects/spatial-genomics")
liver_h5ad <- sc$read_h5ad("./data/all_batches_mouse_only_raw_proximities.h5ad")
liver_h5ad

liver_mat_raw <- liver_h5ad$X
dim(liver_mat_raw) # 5434 cells by 55287 genes (human + mouse genes)
liver_cell_names <- liver_h5ad$obs_names$to_list()
rownames(liver_mat_raw) <- liver_cell_names <- make.names(liver_cell_names, allow_ = F)
colnames(liver_mat_raw) <- liver_gene_names <- liver_h5ad$var_names$to_list()
liver_obs <- liver_h5ad$obs
liver_obs$cell_iden <- rownames(liver_obs)
table(liver_obs$CellType)
liver_mat <- liver_mat_raw[, grep("mm10", colnames(liver_mat_raw))]
# filter out human transcripts
dim(liver_mat) # 5434 cells, 24917 genes

# load location info 
liver_obs$location <- paste0(liver_obs$batch, "_", liver_obs$X, "_", liver_obs$Y)
table(liver_obs$CellType)

# sum(rownames(liver_mat) == liver_obs$cell_iden)
liver_tumor_mat <- liver_mat[liver_obs$CellType == "mc38", ]
dim(liver_tumor_mat) # 2256 mc38 cells and 24917 genes

# normalization
# multiplied each entry by (median / current_sum) so that the sum across the row is equal to the median.
med_sum <- median(liver_obs$n_counts)
# med_sum <- median(rowSums(liver_tumor_mat)) #median of the sum of counts across the cells
mat_normalized_t <- apply(t(liver_tumor_mat), 2, function(i) i * med_sum / sum(i))
liver_tumor_mat_norm <- t(mat_normalized_t)
# sum(rowSums(mat_normalized) == med_sum)
dim(liver_tumor_mat_norm) # 2256 cells by 24917 genes (normalized across all the cells)


# 2. identify tumor-specific gene modules using NMF 
# - center the expression matrix individually by removing the mean expression for each gene.
# - set negative values to zero.
# - perform sparse nonsmooth NMF using nsNMF function in NMF package (rank = 20). 

# liver-NMF
liver_tumor_mat_center <- liver_tumor_mat_norm - colMeans(liver_tumor_mat_norm)[col(liver_tumor_mat_norm)]
# range(colMeans(liver_tumor_mat_center))
liver_tumor_mat_nn <- liver_tumor_mat_center
liver_tumor_mat_nn[liver_tumor_mat_nn < 0] <- 0
dim(liver_tumor_mat_nn) # 2256 cells and 24917 genes
liver_tumor_mat_nn_sel <- liver_tumor_mat_nn[, colSums(liver_tumor_mat_nn) != 0]
dim(liver_tumor_mat_nn_sel) #  2256 cells by 20787 genes

# TEST CASE: 
# sample_raw <- liver_tumor_mat_nn_sel[1:1000, 1:2000]
# sample <- sample_raw[rowSums(sample_raw) != 0, colSums(sample_raw) != 0]
# dim(sample)
# # estim.r <- nmf(x = t(sample),
# #                rank = 2:5, nrun=11, seed=123456)
# # png(file = "./estim.png", width = 1024, height = 768)
# # plot(estim.r)
# # dev.off()
# fit1 <- estim.r$fit$`2`
# s_max_fit1 <- extractFeatures(fit1, method = "max")
# s_kim_fit1 <- extractFeatures(fit1)
# lapply(s_max_fit1, write, "./liver_module_max_fit1.txt", append = TRUE, ncolumns=1000)
# lapply(s_kim_fit1, write, "./liver_module_kim_fit1.txt", append = TRUE, ncolumns=1000)
# fit2 <- nmf(x = t(sample), # feature by sample
#             rank = 10, method = "nsNMF")
# s_max_fit2 <- extractFeatures(fit2, method = "max")
# s_kim_fit2 <- extractFeatures(fit2)
# lapply(s_max_fit2, write, "./liver_module_max_fit2.txt", append = TRUE, ncolumns=1000)
# lapply(s_kim_fit2, write, "./liver_module_kim_fit2.txt", append = TRUE, ncolumns=1000)

nmf_res_liver <- nmf(x = t(liver_tumor_mat_nn_sel), # feature by sample
                    rank = 10, method = "nsNMF")
nmf_res_liver@fit
# <Object of class:NMFns>
# features: 764 
# basis/rank: 10
# samples: 2256
# theta: 0.5 
s_max <- extractFeatures(nmf_res_liver, method = "max")
s_kim <- extractFeatures(nmf_res_liver)

lapply(s_max, write, "./liver_module_max.txt", append = TRUE, ncolumns=5000)
lapply(s_max, cat, "\n", file = "./liver_module_max.txt", append = TRUE)
# lapply(s_kim, write, "./liver_module_kim.txt", append = TRUE, ncolumns=5000)
