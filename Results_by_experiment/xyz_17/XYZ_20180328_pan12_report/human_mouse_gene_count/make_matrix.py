import csv


def make_dict(map):

    map_dict = {}

    with open(map) as f:
        for line in f:
            barcode, x, y = line.split("\t")
            map_dict[barcode] = (x, y)
    return map_dict


def annotate_xy(barcode_annotate, map_dict):

    annotate_dict = {}

    with open(barcode_annotate) as f:
        for line in f:
            id, num = line.split(",")
            pcr, bar = id.split(".")
            
            annotate_dict[num.strip("\n")] = (map_dict[bar], pcr)

    return annotate_dict

def annotate_counts(counts, annotated_dict):

    counts_annotated = {}

    with open(counts) as f:
        for line in f:
            gene_id, num, count = line.split(",")
            
            counts_annotated[gene_id, num] = (annotated_dict[num], count.strip("\n"))

    return counts_annotated

if __name__ == "__main__":

    map="../plate3_barcode_list.txt"

    barcode_num = "cell_annotate.txt"

    map_dict = make_dict(map)
    annotated_dict = annotate_xy(barcode_num, map_dict)

    counts = "count.MM"

    output = annotate_counts(counts, annotated_dict)

    w = csv.writer(open("spatial_output.csv", "w"))
    for key, val in output.items():
        w.writerow([key[0], key[1], val[0][0][0].strip("\n"), val[0][0][1].strip("\n"), val[0][1], val[1]])
