from collections import defaultdict
import pickle
import csv

cells_dict = defaultdict(dict)

genes = set()

with open("count.MM") as f:
	for line in f:
		gene, cell, value = line.split(",")
		genes.add(gene)
		cells_dict[cell][gene] = value

pickle.dump(cells_dict, open("cells_dict.pickle", "wb"))

genes = sorted(list(genes))

with open("output.csv", "w") as out:
	out_csv = csv.writer(out, delimiter="\t")
	csv.writerow(["gene:"]+genes)
	for cell_id, cell_dict in cells_dict:
		csv.writerow([cell_id]+[cell_dict.get(gene, 0) for gene in genes])
