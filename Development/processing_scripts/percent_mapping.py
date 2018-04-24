# find the % mapping for all samples in directory provided 
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

import glob
import matplotlib as plt
import csv


dir = "/ye/yelabstore2/spatial/XYZ_20180122_456_2/STAR_alignment/"
files = glob.glob(dir + "*Log.final.out")

percent_mapping = []


for file in files:
    with open(file) as f:
        for line in f:
            split_line = line.split("\t")
            if " Uniquely mapped reads %" in split_line[0]:
                percent_mapping.append(split_line[1].strip("\n").strip("%"))




plt.barplot(percent_mapping, col="violet")
plt.save("mapping_percentages.png")

report_name = dir + "percent_mapping"
#
# with open(report_name, "w") as f:
#     for item in per