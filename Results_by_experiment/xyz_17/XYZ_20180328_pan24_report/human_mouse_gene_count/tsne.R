# Using seurat to perform tSNE unsupervised clustering 
# XYZeq 
# Ye and Marsen labs @ UCSF 
# March 2018 
# author <christa.caggiano@ucsf.edu> 

# requires packages Rtsne, ggplot2, Seurat, and irlba 

# install.packages("Rtsne")
# install.packages("ggplot2")
# install.packages("Seurat")
# install.packages("irlba")

# clear memory 
rm(list=ls())


library("Rtsne")
library("ggplot2")
library("Seurat")
library("irlba")



  ##########################################################################

gene_count = read.csv("spatial_output_concat.csv", header=T)

row.names(gene_count)=gene_count$X
gene_count = subset(gene_count, select = -c(X) )
min(colSums (gene_count, na.rm = FALSE, dims = 1))
sd(colSums (gene_count, na.rm = FALSE, dims = 1))

pan = CreateSeuratObject(raw.data = gene_count, min.cells = 3, project="pancreas") 
                           
par(mfrow = c(1, 1))

VlnPlot(object = pan, features.plot = c("nGene", "nUMI"), nCol = 2)

pan <- FilterCells(object = pan, subset.names = c("nGene"), low.thresholds = c(200), high.thresholds = c(2500))
GenePlot(object = pan, gene1 = "nUMI", gene2 = "nGene")

pan <- NormalizeData(object = pan, normalization.method = "LogNormalize", 
                      scale.factor = 10000)
pan <- FindVariableGenes(object = pan, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pan@var.genes)
pan <- ScaleData(object = pan, vars.to.regress = c("nUMI"))

pan <- RunPCA(object = pan, pc.genes = pan@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)

VizPCA(object = pan, pcs.use = 1:2)

PCAPlot(object = pan, dim.1 = 1, dim.2 = 2)

PCHeatmap(object = pan, cells.use = 200, pc.use = 1, do.balanced = TRUE, label.columns = FALSE)

PCHeatmap(object = pan, pc.use = 1:12, cells.use = 200, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

pan <- JackStraw(object = pan, num.replicate = 100, display.progress = FALSE)

JackStrawPlot(object = pan, PCs = 1:12)
PCElbowPlot(object = pan)

pbmc <- FindClusters(object = pan, reduction.type = "pca", dims.use = 1:20, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = pan)


pan <- RunTSNE(object = pan, dims.use = 1:20, do.fast = TRUE)
TSNEPlot(object = pan)

cluster1.markers <- FindAllMarkers(object = pan)









