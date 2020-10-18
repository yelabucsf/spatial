python_dir <- "/Applications/anaconda3/bin/python3"
working_dir <- "~/Documents/projects/spatial-genomics"
setwd(working_dir) # set your working directory
# your working directory should contain a `./data`` folder and a `./result` folder to save figures and csvs.
result_dir <- "./result"
data_dir <- "./data"

# 1. load packages
library("reticulate")
use_python(python_dir, required = T)
sc <- import("scanpy")
library(tidyverse)
library(NMF)
install.extras('NMF')

# 2. load liver data (decontaminated and preprocessed - selecting cells)
decont_mat <- read.csv(paste0(data_dir, "/decont_mat_after_v2.csv"), 
                       row.names = 1)
cell_names_raw <- make.names(colnames(decont_mat), allow_ = F)

liver_h5ad <- sc$read_h5ad(paste0(data_dir, "/all_batches_mouse_only_raw_proximities.h5ad"))
liver_h5ad # 5434 cells by 55287 genes
liver_cell_names <- make.names(liver_h5ad$obs_names$to_list(), allow_ = F)
liver_obs <- liver_h5ad$obs
liver_obs$cell_iden <- rownames(liver_obs)
# table(liver_obs$CellType)

liver_mat <- decont_mat[grep("mm10", rownames(decont_mat)), # mouse genes only
                        colnames(decont_mat) %in% liver_cell_names] # selected cells only
dim(liver_mat) # 24917 genes by 5434 cells

# 3. filter to retain genes with high expression
genes_to_keep <- which(rowSums(edgeR::cpm(liver_mat) > 5) >= 15)
liver_mat_filter <- liver_mat[genes_to_keep, ]
dim(liver_mat_filter) # 14222 genes by 5434 cells

# 4. normalization
med_sum <- median(colSums(liver_mat_filter)) #median of the sum of counts across the cells = 377.4436
mat_normalized_t <- apply(liver_mat_filter, 2, function(i) i * med_sum / sum(i))
liver_mat_norm <- t(mat_normalized_t)
# sum(rowSums(liver_mat_norm) == med_sum)
dim(liver_mat_norm) # 5434 cells by 14222 genes (normalized across all the cells)

# 5. select mc38 cells
liver_tumor_mat_norm <- liver_mat_norm[liver_obs$CellType == "mc38", ]
dim(liver_tumor_mat_norm) # 2256 mc38 cells and 14222 genes

# 6. identify tumor-specific gene modules using NMF 
# - center the expression matrix individually by removing the mean expression for each gene.
# - set negative values to zero.
# - perform sparse nonsmooth NMF using nsNMF function in NMF package (rank = 20). 

# liver-NMF
liver_tumor_mat_center <- liver_tumor_mat_norm - colMeans(liver_tumor_mat_norm)[col(liver_tumor_mat_norm)]
# range(colMeans(liver_tumor_mat_center))
liver_tumor_mat_nn <- liver_tumor_mat_center
liver_tumor_mat_nn[liver_tumor_mat_nn < 0] <- 0
dim(liver_tumor_mat_nn) # 2256 cells and 14222 genes
liver_tumor_mat_nn_sel <- liver_tumor_mat_nn[, colSums(liver_tumor_mat_nn) != 0]
dim(liver_tumor_mat_nn_sel) #  2256 cells by 14192 genes

estim.r <- nmf(x = t(liver_tumor_mat_nn_sel),
               rank = 10:20, nrun = 20, method = "nsNMF")
png(file = paste0(result_dir, "/estim.png"), width = 1024, height = 768)
plot(estim.r)
dev.off()

for(i in 10:20){
  fit_i <- estim.r$fit[[i-9]]
  # fit1 <- estim.r$fit$`i`
  s_max_fit_i <- extractFeatures(fit_i, method = "max")
  lapply(s_max_fit_i, write, paste0(result_dir, "/liver_module_max_write", i, ".txt"), append = TRUE, ncolumns=1000)
  lapply(s_max_fit_i, cat, "\n", file = paste0(result_dir, "/liver_module_max_cat", i, ".txt"), append = TRUE)
}


# TEST CASE: 
# sample_raw <- liver_tumor_mat_nn_sel[1:100, 1:2000]
# sample <- sample_raw[rowSums(sample_raw) != 0, colSums(sample_raw) != 0]
# dim(sample)
# estim.r <- nmf(x = t(sample),
#                rank = 2:5, nrun=11, seed=123456)
# png(file = paste0(result_dir, "/estim.png"), width = 1024, height = 768)
# plot(estim.r)
# dev.off()
# for(i in 2:5){
#   fit_i <- estim.r$fit[[i-1]]
#   # fit1 <- estim.r$fit$`i`
#   s_max_fit_i <- extractFeatures(fit_i, method = "max")
#   lapply(s_max_fit_i, cat, "\n", file = paste0(result_dir, "/liver_module_max_fit", i, ".txt"), append = TRUE)
#   lapply(s_max_fit_i, write, paste0(result_dir, "/liver_module_max_write", i, ".txt"), append = TRUE, ncolumns=1000)
# }
# s_kim_fit1 <- extractFeatures(fit1)
# lapply(s_max_fit1, write, "./liver_module_max_fit1.txt", append = TRUE, ncolumns=1000)
# lapply(s_kim_fit1, write, "./liver_module_kim_fit1.txt", append = TRUE, ncolumns=1000)
# fit2 <- nmf(x = t(sample), # feature by sample
#             rank = 20, method = "nsNMF")
# s_max_fit2 <- extractFeatures(fit2, method = "max")
# s_kim_fit2 <- extractFeatures(fit2)
# lapply(s_max_fit2, write, "./liver_module_max_fit2.txt", append = TRUE, ncolumns=1000)
# lapply(s_kim_fit2, write, "./liver_module_kim_fit2.txt", append = TRUE, ncolumns=1000)
# lapply(s_max_fit2, cat, "\n", file = "./liver_module_max.txt", append = TRUE)
