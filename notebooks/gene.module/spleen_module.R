# load packages
library("reticulate")
use_python("/Applications/anaconda3/bin/python3", required = T)
sc <- import("scanpy")
library(tidyverse)
library(NMF)
install.extras('NMF')

# load spleen data (without decontamination yet preprocessed)
setwd("~/Documents/Github/spatial-genomics")
spleen_h5ad <- sc$read_h5ad("./data/spleen.concat.raw_counts.h5ad")
spleen_h5ad
# AnnData object with n_obs × n_vars = 7505 × 52025 
# obs: 'batch', 'percent_mito_mouse', 'percent_mito_human', 'n_counts', 'n_genes', 'mouse_counts', 'human_counts', 'ratio', 'cell_call', 'leiden', 'CellType'
# var: 'n_cells', 'log_cells'

spleen_mat_raw <- spleen_h5ad$X
dim(spleen_mat_raw) # 7505 cells by 52025 genes (human + mouse genes)
spleen_cell_names <- spleen_h5ad$obs_names$to_list()
rownames(spleen_mat_raw) <- spleen_cell_names <- make.names(spleen_cell_names, allow_ = F)
colnames(spleen_mat_raw) <- spleen_gene_names <- spleen_h5ad$var_names$to_list()
spleen_obs <- spleen_h5ad$obs
spleen_obs$cell_iden <- make.names(rownames(spleen_obs), allow_ = F)
table(spleen_obs$CellType)
spleen_mat <- spleen_mat_raw[, grep("mm10", colnames(spleen_mat_raw))]
# filter out human transcripts
dim(spleen_mat) # 7505 cells, 22810 genes

# load location info 


# sum(rownames(spleen_mat) == spleen_obs$cell_iden) # 7505
spleen_tumor_mat <- spleen_mat[spleen_obs$CellType == "mc38", ]
dim(spleen_tumor_mat) # 2356 mc38 cells and 22810 genes

# normalization
# multiplied each entry by (median / current_sum) so that the sum across the row is equal to the median.
med_sum <- median(spleen_obs$n_counts)
# med_sum <- median(rowSums(spleen_tumor_mat)) #median of the sum of counts across the cells
mat_normalized_t <- apply(t(spleen_tumor_mat), 2, function(i) i * med_sum / sum(i))
spleen_tumor_mat_norm <- t(mat_normalized_t)
# sum(colSums(mat_normalized_t) == med_sum) # 2356
dim(spleen_tumor_mat_norm) # 2356 cells by 22810 genes (normalized across all the cells)


# 2. identify tumor-specific gene modules using NMF 
# - center the expression matrix individually by removing the mean expression for each gene.
# - set negative values to zero.
# - perform sparse nonsmooth NMF using nsNMF function in NMF package (rank = 20). 

# spleen-NMF
spleen_tumor_mat_center <- spleen_tumor_mat_norm - colMeans(spleen_tumor_mat_norm)[col(spleen_tumor_mat_norm)]
# range(colMeans(spleen_tumor_mat_center))
spleen_tumor_mat_nn <- spleen_tumor_mat_center
spleen_tumor_mat_nn[spleen_tumor_mat_nn < 0] <- 0
dim(spleen_tumor_mat_nn) # 2356 cells and 22810 genes
spleen_tumor_mat_nn_sel <- spleen_tumor_mat_nn[, colSums(spleen_tumor_mat_nn) != 0]
dim(spleen_tumor_mat_nn_sel) #  2356 cells by 18898 genes

# TEST CASE: 
# sample_raw <- spleen_tumor_mat_nn_sel[1:100, 1:200]
# sample <- sample_raw[rowSums(sample_raw) != 0, colSums(sample_raw) != 0]
# dim(sample)
# # estim.r <- nmf(x = t(sample),
# #                rank = 2:5, nrun=11, seed=123456)
# # png(file = "./estim.png", width = 1024, height = 768)
# # plot(estim.r)
# # dev.off()
# # fit1 <- estim.r$fit$`2`
# # s_max_fit1 <- extractFeatures(fit1, method = "max")
# # s_kim_fit1 <- extractFeatures(fit1)
# # lapply(s_max_fit1, write, "./spleen_module_max_fit1.txt", append = TRUE, ncolumns=1000)
# # lapply(s_kim_fit1, write, "./spleen_module_kim_fit1.txt", append = TRUE, ncolumns=1000)
# fit2 <- nmf(x = t(sample), # feature by sample
#             rank = 10, method = "nsNMF")
# s_max_fit2 <- extractFeatures(fit2, method = "max")
# s_kim_fit2 <- extractFeatures(fit2)
# lapply(s_max_fit2, write, "./spleen_module_max_fit2.txt", append = TRUE, ncolumns=1000)
# lapply(s_kim_fit2, write, "./spleen_module_kim_fit2.txt", append = TRUE, ncolumns=1000)

nmf_res_spleen <- nmf(x = t(spleen_tumor_mat_nn_sel), # feature by sample
                     rank = 10, method = "nsNMF")
nmf_res_spleen@fit
# <Object of class:NMFns>
# features: 764 
# basis/rank: 10
# samples: 2256
# theta: 0.5 
s_max <- extractFeatures(nmf_res_spleen, method = "max")
s_kim <- extractFeatures(nmf_res_spleen)

lapply(s_max, write, "./spleen_module_max.txt", append = TRUE, ncolumns=5000)
# lapply(s_kim, write, "./spleen_module_kim.txt", append = TRUE, ncolumns=5000)
