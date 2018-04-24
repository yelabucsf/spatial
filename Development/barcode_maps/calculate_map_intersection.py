# find intersection of two barcode lists, mostly just for sanity checking
# winter 2018 
# author <christa.caggiano@ucsf.edu>  


s1 = set()
with open("primarybarcode34.txt") as f: 
	for line in f: 
		s1.add(line)

s2 = set()
with open("plate34_barcode_list.txt") as f: 
	for line in f: 
		s2.add(line)
		
print("in primary not in secondary")
print(s1-s2)
print("in secondary not in primary")
print(s2-s1)
