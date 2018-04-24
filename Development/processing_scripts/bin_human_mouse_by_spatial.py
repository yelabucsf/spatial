# sum up all the mouse and human UMIs for a given spatial well 
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

# imports 
import csv

# row in annotated output that contains "hg" or "mm" for species 
id_num = 1

# open the annotated UMI counts file 
with open("annotated_output.csv") as f, open("human-mouse-sums.csv", "w") as new_f:
    
    c = csv.reader(f)  # read as a csv 
    d = {}

    for line in c: 
        key = (line[2], line[3]). # set the key as the coordinates of the well 
        if key not in d:
            d[key] = [0,0]

        if line.split()[0] == "hg": # if the line is human, add the human UMI counts 
            d[key][0] += int(line[5])
        else: 
            d[key][1] += int(line[5])  # otherwise, add as mouse 

    # output formatted file 
    out = csv.writer(new_f)
    for key, value in d.items():
        out.writerow(list(key)+value)