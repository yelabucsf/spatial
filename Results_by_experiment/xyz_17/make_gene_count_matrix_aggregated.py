# concatenates several pancreas outputs by spatial well in order to get 
# higher resolution for tsne plotting 
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

# imports 
import csv
from collections import defaultdict
results = defaultdict(dict)

# initializing values to find the biggest well size in particular chip 
max_x = 1
max_y = 1

# open a combined annotated matrix 
with open("pancreas_combined.csv") as f:
    lines = csv.reader(f)

    for line in lines:
        gene_number = line[0] # not used
        
        x = int(line[1]) # get spatial coordinates 
        y = int(line[2])

        value = line[3]  # number of UMIS 

        gene_id = line[4]  # gene number assigned by pipeline 

        results[gene_id][(x,y)] = int(value)  # make a dictionary where spatial coordinates are key and map to UMI value 

        max_x = max(max_x, int(x))  # update the max values, to keep track of chip size 
        max_y = max(max_y, int(y))


# write output 
with open("spatial_output_concat.csv", "w") as f:

    file = csv.writer(f)
    coordinates = []

    # for the max sizes determined above, append a list of coordinates to be filled in 
    for i in range(1, max_x+1):
        for j in range(1, max_y+1):
            coordinates.append((i,j))

    # write the first row of the matrix as the coordinates 
    first_row = [""]+coordinates
    file.writerow(first_row)

    # for the rest of the matrix, write the UMI counts for the specific gene  
    for gene, d in results.items():
        row = [gene]

        for (x, y) in coordinates:

            row.append(d.get((x,y), 0))

        file.writerow(row)
