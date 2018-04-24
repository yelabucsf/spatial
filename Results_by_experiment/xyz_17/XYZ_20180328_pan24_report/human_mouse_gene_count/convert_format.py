import csv
from collections import defaultdict
results = defaultdict(dict)
max_gene_number = 1
with open("spatial_output.csv") as f:
    lines = csv.reader(f)
    for line in lines:
        gene_id = line[-1]
        results[int(line[0])][int(line[1])] = int(line[-1])
        max_gene_number = max(max_gene_number, int(line[1]))

with open("spatial_output_concat.csv", "w") as f:
    file = csv.writer(f)
    first_row = [""]+list(range(1,max_gene_number+1))
    file.writerow(first_row)
    for gene, d in results.items():
        row = [gene]
        for i in range(1, max_gene_number+1):
            row.append(d.get(i, 0))
        file.writerow(row)
