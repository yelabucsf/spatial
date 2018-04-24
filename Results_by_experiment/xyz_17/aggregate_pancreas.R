# Merging Pancreas spatial sequencing data  
# XYZeq 
# Ye and Marsen labs @ UCSF 
# March 2018 
# author <christa.caggiano@ucsf.edu> 

# clear memory 
rm(list = ls())

# set working directory
setwd("~/Documents/UCSF_year1/Ye-rotation2/")

ChangeNames <- function(x) {
  names(x) <- c("gene", "cell", "x", "y", "pcr", "val")
  return(x)
}



# read in pancreas data 
pan12 = read.csv("XYZ_20180328_pan12_report/human_mouse_gene_count/spatial_output.csv", header=F)
pan24 = read.csv("XYZ_20180328_pan24_report/human_mouse_gene_count/spatial_output.csv", header=F)
pan40 = read.csv("XYZ_20180328_pan40_report/human_mouse_gene_count/spatial_output.csv", header=F)


colnames(pan40) = c("gene", "cell", "x", "y", "pcr", "val")


# make a list of all three datasets 
pancreas_list = c(pan12, pan24, pan40)

# standardize the column names of dataframes of interest 
pancreas_list = lapply(pancreas_list, ChangeNames)

gene_annotation = read.csv("XYZ_20180328_pan24_report/human_mouse_gene_count/gene_name_annotate.txt", header = F)
gene_annotation = subset(gene_annotation, gene_annotation$V3 == "exon")
colnames(gene_annotation) = c("name", "type", "cat", "symbol", "gene")

pan12 = merge(pan12, gene_annotation, by="gene")
sum(pan12$val)/length(unique(pan12$cell))

pan24 = merge(pan24, gene_annotation, by="gene")
sum(pan24$val)/length(unique(pan24$cell))

pan40 = merge(pan40, gene_annotation, by="gene")
sum(pan40$val)/length(unique(pan40$cell))

all_pan = merge(pan12, pan24, by=c("gene", "x", "y"), all=T)
all_pan = merge(all_pan, pan40, by=c("gene", "x", "y"), all=T)

all_pan$value = rowSums(all_pan[,c("val.x", "val.y", "val")], na.rm=TRUE)

all_pan$symbol.x <- as.character(all_pan$symbol.x)
all_pan$symbol.x[is.na(all_pan$symbol.x)] = ""

for (row in 1:nrow(all_pan)){
  if (all_pan[row, "symbol.y"] == all_pan[row, "symbol.x"]) {
    all_pan[row, "name"] = all_pan[row, "symbol.x"]
  } else if (all_pan[row, "symbol.y"] == ""){ 
    all_pan[row, "name"] = all_pan[row, "symbol.x"]
  } else {
    all_pan[row, "name"] = all_pan[row, "symbol.y"]
  }
} 

columns_to_keep = c("gene", "x", "y", "value", "name")

matrix = all_pan[, columns_to_keep]
write.csv(matrix, "pancreas_combined.csv", row.names = F)
