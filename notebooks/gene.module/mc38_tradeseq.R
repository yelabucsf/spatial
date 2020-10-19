python_dir <- "/Applications/anaconda3/bin/python3"
working_dir <- "~/Documents/projects/spatial-genomics"
setwd(working_dir) # set your working directory
num_cores <- 8
# 1. load-packages
library("reticulate")
use_python(python_dir, required = T)
sc <- import("scanpy")

library("tidyverse")
library("tradeSeq")
library("zinbwave")
library("BiocParallel")
library("doParallel")

# Liver data is decontaminated by Yutong. 
# Cells filtered by Derek are used in downstream analysis.

# 2. load-data 
# decont liver mat 
# rm(list=setdiff(ls(), c("decont_mat", "decont_mat_raw)))
decont_mat <- read.csv("./data/decont_mat_after_v2.csv", row.names = 1)
# decont_mat <- round(decont_mat_raw)
# # decont_mat <- decont_mat_raw 
dim(decont_mat) # 55287 genes, 10280 cells (the 59 contaminated genes were removed)
cell_names_raw <- make.names(colnames(decont_mat), allow_ = F)
gene_names <- rownames(decont_mat)

# 3. selected cells
liver_h5ad <- sc$read_h5ad("./data/all_batches_mouse_only_raw_proximities.h5ad")
liver_h5ad # 5434 cells by 55287 genes
liver_cell_names <- make.names(liver_h5ad$obs_names$to_list(), allow_ = F)
liver_obs <- liver_h5ad$obs
liver_obs$cell_iden <- rownames(liver_obs)
# table(liver_obs$CellType)

liver_mat <- decont_mat[grep("mm10", rownames(decont_mat)), # mouse genes only
                        colnames(decont_mat) %in% liver_cell_names] # selected cells only
# dim(liver_mat) # 24917 genes by 5434 cells

# 4. filtering to retain genes with high expression:
# among 24917 genes, there are 14222 genes which have at least 5 cells whose 
# counts for that gene is larger than 5.
genes_to_keep <- which(rowSums(edgeR::cpm(liver_mat) > 5) >= 15)
# tradeSeq: gene by cell in cufflinksCountData
liver_mat_filter <- liver_mat[genes_to_keep, ]
# dim(liver_mat_filter) # 14222 genes by 5434 cells

# 5. normalization
med_sum <- median(colSums(liver_mat_filter)) #median of the sum of counts across the cells = 860.2559
mat_normalized_t <- apply(liver_mat_filter, 2, function(i) i * med_sum / sum(i))
liver_mat_norm <- t(mat_normalized_t)
# sum(rowSums(liver_mat_norm) == med_sum)
dim(liver_mat_norm) # 5434 cells by 14222 genes (normalized across all the cells)

# 6. select-mc38
# dim(liver_obs) # 5434 cells by 22 annotations
# sum(liver_obs$cell_iden == rownames(liver_mat_norm)) # 5434
liver_mc38_mat_norm <- liver_mat_norm[liver_obs$CellType == "mc38", ]
# dim(liver_mc38_mat_norm) # 2256 cells by 14222 genes
liver_mc38_obs <- liver_obs[liver_obs$CellType == "mc38", ]
# table(liver_mc38_obs$batch)
liver_mc38_mat_unnorm <- round(t(liver_mat_filter)[liver_obs$CellType == "mc38", ])
# dim(liver_mc38_mat_unnorm) # 2256 cells by 14222 genes
# range(liver_mc38_mat_unnorm) # 0-224

# which gene has zero counts in mc38
liver_mc38_mat_unnorm_filter <- liver_mc38_mat_unnorm[, which(colSums(liver_mc38_mat_unnorm) != 0)]
dim(liver_mc38_mat_unnorm_filter) # 2256 cells by 14181 genes

# 7. ZI-weights-zinbwave
registerDoParallel(num_cores)
register(DoparParam())
# table(liver_obs$bins_1)
sum_exp_obj <- SummarizedExperiment(t(liver_mc38_mat_unnorm_filter), 
                                    colData = data.frame(bin = liver_mc38_obs$bins_1, 
                                                         batch = liver_mc38_obs$batch))
zinb_res <- zinbFit(sum_exp_obj, X = '~ bin + batch', 
                    commondispersion = TRUE)

save(zinb_res, file = "./mc38_zinb_fit.rda")
zinb_weights <- computeObservationalWeights(zinb_res, t(liver_mc38_mat_unnorm_filter))

# 8. tradeSeq
gam_list <- tradeSeq::fitGAM(t(liver_mc38_mat_unnorm_filter),
                             U = model.matrix(~ - 1 + liver_mc38_obs$batch + liver_mc38_obs$bins_1), 
                             pseudotime = liver_mc38_obs$prox_2, 
                             cellWeights = rep(1, nrow(liver_mc38_obs)), 
                             weights = zinb_weights, 
                             nknots = 6)
save(gam_list, file = "./mc38_tradeseq_list.rda")