python_dir <- "/Applications/anaconda3/bin/python3"
working_dir <- "~/Documents/projects/spatial-genomics"
setwd(working_dir) # set your working directory
num_cores <- 3
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

# decont_mat <- read.csv("./data/decont_mat_after_v2.csv", row.names = 1)
# save(decont_mat, file = "./data/decont_mat.rda")
load("./data/decont_mat.rda")

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
liver_mat_filter_raw <- liver_mat[genes_to_keep, ]
# dim(liver_mat_filter) # 14222 genes by 5434 cells

a <- liver_mat_filter_raw[rownames(liver_mat_filter), liver_obs$CellType == "mc38"]
dim(a) # 2000 genes and 2256 cells
sum(colSums(a) == 0)


# 4.5 select top 2000 most variable genes
# convert to log 2 counts per million
liver_mat_log2 <- edgeR::cpm(liver_mat_filter_raw, log = TRUE)
dim(liver_mat_log2) # 14222 by 5434
Calc_CV <- function(x){sd(x) / mean(x)}
liver_gene_cv <- apply(liver_mat_log2, 1, Calc_CV)
summary(liver_gene_cv)
liver_mat_filter <- liver_mat_log2[which(rank(liver_gene_cv) > length(liver_gene_cv) - 2000), ]
dim(liver_mat_filter) # 2000 genes by 5434 cells

write.csv(rownames(liver_mat_filter), file = "./data/top2000genes.csv")
