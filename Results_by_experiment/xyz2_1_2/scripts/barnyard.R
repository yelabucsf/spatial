# making barnyard plots for cell tiling experiments 
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

# clear memory and set working directory 
rm(list=ls())
setwd("~/Desktop/Ye_lab_desktop/")

# packages 
library('ggplot2')

#############################  plot barnyard read counts ######################################

# read in reads counts for different plate orientations 
counts = read.csv("reads_counts_all.csv", header=T)
spatial = read.csv("reads_spatial.csv", header=F)
spatial_x = read.csv("reads_spatial_x.csv", header=F)
spatial_y = read.csv("reads_spatial_y.csv", header=F)
spatial_pcr = read.csv("reads_spatial_pcr.csv", header=F)

# plot counts for single cells (spatial + PCR)
ggplot(counts,aes(x=human,y=mouse)) +
  geom_point( size=2, color="aquamarine4") + 
  labs(x = "Single Cell", y="Counts") + 
  theme_light() + 
  scale_colour_brewer(palette = "BrBG") + 
  theme(text = element_text(family = "Tahoma", color = "grey20"))

# plot spatial barcode reads counts
ggplot(spatial,aes(x=V3,y=V4)) +
  geom_point( size=2, color="aquamarine3") + 
  labs(x = "Spatial", y="Counts") + 
  theme_light() + 
  scale_colour_brewer(palette = "BrBG") + 
  theme(text = element_text(family = "Tahoma", color = "grey20"))

# plot counts by x-direction (columns) of chip 
ggplot(spatial_x,aes(x=V2,y=V3)) +
  geom_point( size=2, color="aquamarine2") + 
  labs(x = "X-columns", y="Counts") + 
  theme_light() + 
  scale_colour_brewer(palette = "BrBG") + 
  theme(text = element_text(family = "Tahoma", color = "grey20"))

# plot counts in y-direction (rows) of chip 
ggplot(spatial_y,aes(x=V2,y=V3)) +
  geom_point( size=2, color="aquamarine1") + 
  labs(x = "Y-columns", y="Counts") + 
  theme_light() + 
  scale_colour_brewer(palette = "BrBG") + 
  theme(text = element_text(family = "Tahoma", color = "grey20"))

# plot counts only by PCR barcode 
ggplot(spatial_pcr,aes(x=V2,y=V3)) +
  geom_point( size=2, color="darkseagreen") + 
  labs(x = "PCR well", y="Counts") + 
  theme_light() + 
  scale_colour_brewer(palette = "BrBG") + 
  theme(text = element_text(family = "Tahoma", color = "grey20"))


#############################  UMI distribution by species ######################################

# read in UMI counts with human/mouse annotation
ann = read.csv("annotated_output.csv", header=F)

# subset to human and mouse dfs 
ann_hum = subset(ann, ann$V1=="hg")
ann_mouse = subset(ann, ann$V1=="mm")

# calculate number of UMIs in each species
sum(ann_mouse$V4)
sum(ann_hum$V4)

# number of cells for each species 
length(unique(ann_mouse$V3))
length(unique(ann_hum$V3))

# aggregate average number of UMIs by cell for mouse and human  
agg_mouse = aggregate(ann_mouse$V4,list(ann_mouse$V3),mean)
agg_hum = aggregate(ann_hum$V4,list(ann_hum$V3),mean)

# number of mouse UMIs vs human UMIs
barplot(agg_mouse$x, col="red")
barplot(table(agg_hum$x), add = T, col="blue")

# distribution of the human and mouse UMI counts
boxplot(agg_mouse$x, agg_hum$x)

#############################  UMI distribution by species ######################################

# read in spatial output w/ UMI counts by species 
sp = read.csv("spatial_output.csv", header=F)

# subset to human/mouse 
sp_human = subset(sp, sp$V2=="hg")
sp_mouse = subset(sp, sp$V2=="mm")

# aggregate UMI counts by cell number  
sp_mouse_agg = aggregate(sp_mouse$V7,list(sp_mouse$V3),sum)
sp_hum_agg = aggregate(sp_human$V7,list(sp_human$V3),sum)

# make into one dataframe because ggplot is so so dumb 
merged = merge(sp_mouse_agg, sp_hum_agg, by.x=c("Group.1"), by.y=c("Group.1"), all=T)

# plot mouse UMIs and human UMIs by single cells (PCR + spatial barcodes)
ggplot(merged,aes(x=x.y,y=x.x)) +
  geom_point( size=2, color="aquamarine4") + 
  theme_light() + 
  scale_colour_brewer(palette = "BrBG") + 
  theme(text = element_text(family = "Tahoma", color = "grey20"))











