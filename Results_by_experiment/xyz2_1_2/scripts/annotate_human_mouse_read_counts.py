# takes human - mouse read count fraction and maps onto spatial x/y coordinates 
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

#imports 
import csv

def make_dict(map):
    """
    @param map: tab separated file containing spatial coordinates 
    mapping to identifying barcodes

    return map_dict: dictionary of {barcode: (x,y)} 
    """

    map_dict = {}

    with open(map) as f:

        for line in f:
            barcode, x, y = line.split("\t")
            map_dict[barcode] = (x, y) # for each barcode, store its x/y coordinates

    return map_dict


def annotate_xy(barcode_annotate, map_dict):
    """
    @param map_dict: dictionary containing barcodes and their xy coordinates
    @param barcode_annotate: file to be annotated, in this case, human/mouse read counts 

    return annotate_dict: {PCR.barcode:(x, y), human_counts, mouse_counts, celegan_counts}  
    """

    annotate_dict = {}

    with open(barcode_annotate) as f:

        for line in f:

            id, human, mouse, celegan = line.split(",")  # split line by columns (celegan should be 0 for our experiments)
            pcr, bar = id.split(".") # split identifier into its components, spatial barcode and PCR barcode
            annotate_dict[id] = (map_dict[bar], human, mouse, celegan) # find the coordinates for the given barcode 


    return annotate_dict


if __name__ == "__main__":

    # not a great way of doing this, but read in this hardcoded barcode map file 
    map="plate34_barcode_list.txt"
    
    # make the dictionary 
    map_dict = make_dict(map)

    # read in the human mouse read counts 
    human_mouse = "human_mouse_fraction.txt"

    # annotate the human mouse barcode list 
    annotated_dict = annotate_xy(human_mouse, map_dict)

    # output as a csv 
    w = csv.writer(open("human-mouse.csv", "w"))
    for key, val in annotated_dict.items():
        w.writerow([key, val[0][0], val[0][1].strip("\n"), val[1], val[2]])
