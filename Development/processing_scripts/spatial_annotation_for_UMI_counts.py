# generalized script to annotate UMI count information with its spatial position 
# given a barcode map 
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

# imports 
import csv


def make_dict(map):
    """
    @param map: text file containing the barcode sequence and its x, y coordinate on the chip 
    
    return map_dict: {barcode, (x, y)}

    """

    map_dict = {}

    with open(map) as f:

        for line in f:  # each line corresponds to one barcode 

            barcode, x, y = line.split("\t")  # split line into barcode and its x y coords
            map_dict[barcode] = (x, y)

    return map_dict


def annotate_xy(barcode_annotate, map_dict):
    """
    extra step of matching barcode information to the cell id # assigned in the pipeline 

    @param barcode_annotate: file that matches cell ID # to a barcode 
    @ param map_dict: {barcode:(x, y)} 

    return annotate_dict: {cell_number : (x, y), PCR}
    """

    annotate_dict = {}

    with open(barcode_annotate) as f:
        for line in f:

            id, num = line.split(","). # split into barcode and cell ID 
            pcr, bar = id.split(".")  # split into spatial barcode and PCR barcode 
            annotate_dict[num.strip("\n")] = (map_dict[bar], pcr)

    return annotate_dict

def annotate_counts(counts, annotated_dict):

    """
    for each UMI, assign spatial information 
    
    @param counts: UMI counts 
    @param annotated_dict: {cell_number : (x, y), PCR}

    return counts_annotated: {gene: ((x, y), pcr), counts}

    """


    counts_annotated = {}

    with open(counts) as f:

        for line in f:

            gene_id, num, count = line.split(","). # gene id #, cell #, UMI counts
            counts_annotated[gene_id] = (annotated_dict[num], count.strip("\n"))

    return counts_annotated

if __name__ == "__main__":

    # hardcoded barcode file location containing spatial barcode and x/y coordinates  
    map="/ye/yelabstore2/spatial/plate3_barcode_list.txt"

    # file containing cell ID and assigns barcodes (Barcode, ID) 
    barcode_num = "cell_annotate.txt"

    # create map of barcodes:(x,y) coordinates
    map_dict = make_dict(map)

    # annotate the cell ID with the spatial coordinate information 
    annotated_dict = annotate_xy(barcode_num, map_dict)

    # file containing UMI count information 
    counts = "count.MM"

    # annotate the counts 
    output = annotate_counts(counts, annotated_dict)

    # output as CSV 
    w = csv.writer(open("spatial_output.csv", "w"))
    for key, val in output.items():
        w.writerow([key, val[0][0][0], val[0][0][1], val[0][1], val[1]])
