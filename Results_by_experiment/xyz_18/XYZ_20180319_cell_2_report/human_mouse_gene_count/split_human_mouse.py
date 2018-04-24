from __future__ import print_function
import csv

gene_dict = {}

with open("gene_name_annotate.txt") as g: 
	for line in g: 
		gene_name = line.split(",")[0]
		gene_number = line.split(",")[4].strip("\n")
		gene_dict[gene_number] = gene_name


a = open("annotated_output.csv", "w")
with open("count.MM") as s: 
	for line in s: 
		if "hg" in gene_dict[line.split(",")[0]]:
			print("hg," + line.strip("\n"), file=a)
		else: 
			print("mm," +line.strip("\n"), file=a)

