import csv


def make_dict(map):

    map_dict = {}

    with open(map) as f:
        for line in f:
            barcode, x, y = line.split(",")
            map_dict[barcode] = (x, y)
    return map_dict


def annotate_xy(barcode_annotate, map_dict):

    annotate_dict = {}

    with open(barcode_annotate) as f:
        for line in f:
            id, num = line.split(",")
            pcr, bar = id.split(".")
            try: 
                annotate_dict[num.strip("\n")] = (map_dict[bar], pcr)
            except KeyError: 
                pass 
    return annotate_dict
    

def annotate_counts(counts, map_dict):

    counts_annotated = {}

    with open(counts) as f:
        for line in f:
            bar, human, mouse = line.split(",")
            pcr, spatial = bar.split(".")
            try: 
                counts_annotated[mouse.strip("\n"), human.strip("\n")] = (map_dict[spatial.strip("\n")], pcr)
            except KeyError: 
                pass 
    return counts_annotated


if __name__ == "__main__":

    map="../plate23_barcode_map.txt"

    counts = "human_mouse_fraction.txt"

    map_dict = make_dict(map)
    # annotated_dict = annotate_xy(counts, map_dict)

    # # counts = "annotated_output.csv"

    output = annotate_counts(counts, map_dict)

    w = csv.writer(open("spatial_output.csv", "w"))
    for key, val in output.items():
        w.writerow([val[0][0].strip("\n"), val[0][1].strip("\n"), val[1].strip("\n"), key[0].strip("\n"), key[1].strip("\n")])
