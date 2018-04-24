import csv
from collections import defaultdict
results = defaultdict(dict)
max_x = 1
max_y = 1
with open("spatial_output.csv") as f:
    lines = csv.reader(f)
    for line in lines:
        gene_number = line[0] # not used
        x = line[1]
        y = line[2]
        value = line[3]
        gene_id = line[4]
        results[gene_id][(x,y)] = int(value)
        max_x = max(max_x, int(x))
        max_y = max(max_y, int(y))

with open("spatial_output_concat.csv", "w") as f:
    file = csv.writer(f)
    coordinates = []
    for i in range(1, max_x+1):
        for j in range(1, max_y+1):
            coordinates.append((i,j))
    first_row = [""]+list(range(1,max_gene_number+1))
    file.writerow(first_row)
    for gene, d in results.items():
        row = [gene]
        for (x, y) in coordinates
            row.append(d.get((x,y), 0))
        file.writerow(row)
