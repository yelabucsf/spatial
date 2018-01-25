# plotting mapping % 
# jan 2018 
# author <christa.caggiano@ucsf.edu> 

import os

def find_percent_mapped(log_file):
	if os.path.isfile(log_file):
		with open(log_file) as f: 
			for line in f: 
				if "|" in line:
					#no_newline = line.strip(" ")
					title, value = line.split(" |\t")
					if title.strip().startswith("% of reads unmapped: too short"): 
						return value 
	else: 
		return 0 

if __name__ == "__main__": 
	# phix_percent = []
	mm10_percent = [] 

	with open("file_list2.txt") as f: 
		for file_dir in f: 
			aligned_dir = file_dir.strip("\n")
			print(aligned_dir)
			mm10 = aligned_dir + "/Log.final.out"
			# phix = aligned_dir + "Log.final.phix"
			mm10_percent.append(find_percent_mapped(mm10))
			# phix_percent.append(find_percent_mapped(phix))

	with open('percents.txt', 'w') as percent_list:
		for percent in mm10_percent:
			percent_list.write("%s " % percent)
		# percent_list.write('phix')
		# for percent in phix_percent: 
		# 	percent_list.write("%s " % percent)



