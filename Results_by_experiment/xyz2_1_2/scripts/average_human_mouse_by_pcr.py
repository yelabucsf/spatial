# for each PCR barcode, report the average number of human and mouse reads
# winter 2018 
# author <christa.caggiano@ucsf.edu> 

# imports 
import csv
from collections import defaultdict


def avg(l): 
	"""
	calculates average of a list 
	"""
	return sum(l)/len(l)


if __name__=="__main__": 

	# open file containing human mouse read counts for a PCR + spatial barcode 
	with open("human_mouse_fraction.txt") as f, open("result.csv", "w") as out:
		
		d_human = defaultdict(list)
		d_mouse = defaultdict(list)
		d_cele = defaultdict(list)
		ds = [d_human,d_mouse]
		csv_keys = ["human_reads","mouse_reads"]

		# Makes the file into a list of dictionaries
		lines = csv.DictReader(f, delimiter="\t")
		for line in lines:
			key = line['sample'].split(".")[0]
			for csv_key, d in zip(csv_keys,ds):
				d[key] += [int(line[csv_key])]

		# print average results by PCR barcode to csv 
		fieldnames = ['PCR Barcode', 'Human Average', 'Mouse Average', "Cele Average"]
		writer = csv.DictWriter(out, fieldnames=fieldnames)
		for key in d_human:
			writer.writerow(
				{
				'PCR Barcode': key,
				'Human Average': avg(d_human[key]),
				'Mouse Average': avg(d_mouse[key])
				})


