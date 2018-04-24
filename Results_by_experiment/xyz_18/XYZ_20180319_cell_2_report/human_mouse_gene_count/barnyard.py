import csv
if __name__ == "__main__": 


	count_x_y_human = {}
	count_x_y_mouse = {}
	count_y_human = {}
	count_x_human = {}
	count_y_mouse = {}
	count_x_mouse = {}
	count_pcr_mouse = {}
	count_pcr_human = {}
	with open("spatial_output.csv") as f: 
		for line in f: 
			line = line.strip("\n")
			gene, org, cell, x, y, pcr, count = line.split(",")
			x = int(x)
			y = int(y)
			if org == "hg": 
				if x in count_x_y_human: 
					# count_x_y_human[(x,y)] += float(count)
					# count_y_human[x] += float(count)
					count_x_human[y] += float(count)
					# count_pcr_human[pcr] += float(count)


				else: 
					# count_x_y_human[(x,y)] = float(count)
					# count_y_human[x] = float(count)
					count_x_human[y] = float(count)
					# count_pcr_human[pcr] = float(count)

			else: 
				if x in count_x_y_mouse: 
					# count_x_y_mouse[(x,y)] += float(count)
					# count_y_mouse[x] += float(count)
					count_x_mouse[y] += float(count)
					# count_pcr_mouse[pcr] += float(count)
				else: 
					# count_x_y_mouse[(x,y)] = float(count)
					# count_y_mouse[x] = float(count)
					count_x_mouse[y] = float(count)
					# count_pcr_mouse[pcr] = float(count)

		tupule_list = []
		for x in range(1, 18): 
			for y in range(1, 43): 
				tupule_list.append((x,y))

		spatial = {}
		for item in tupule_list: 
			if item in count_x_y_human and item in count_x_y_mouse: 
				spatial[item] = (count_x_y_mouse[item], count_x_y_human[item])
			elif item in count_x_y_human and item not in count_x_y_mouse: 
				spatial[item] = (0, count_x_y_human[item])
			elif item in count_x_y_mouse and item not in count_x_y_human: 
				spatial[item] = (count_x_y_mouse[item], 0)
			else: 
				spatial[item] = (0, 0)


		x_list = [i for i in range(1,18)]

		xx = {}
		for item in x_list: 
			if item in count_x_human and item in count_x_mouse: 
				xx[item] = (count_x_mouse[item], count_x_human[item])
			elif item in count_x_human and item not in count_x_mouse: 
				xx[item] = (0, count_x_human[item])
			elif item in count_x_mouse and item not in count_x_human: 
				xx[item] = (count_x_mouse[item], 0)
			else: 
				xx[item] = (0, 0)

		# y_list = [i for i in range(1,43)]

		# yy = {}
		# for item in x_list: 
		# 	if item in count_y_human and item in count_y_mouse: 
		# 		yy[item] = (count_y_mouse[item], count_y_human[item])
		# 	elif item in count_y_human and item not in count_y_mouse: 
		# 		yy[item] = (0, count_y_human[item])
		# 	elif item in count_y_mouse and item not in count_y_human: 
		# 		yy[item] = (count_y_mouse[item], 0)
		# 	else: 
		# 		yy[item] = (0, 0)

		# pcr = {}
		# for key, val in count_pcr_mouse.items(): 
		# 	pcr[key] = (val, count_pcr_human[key])


	# w = csv.writer(open("umi_xy.csv", "w"))
	# for key, val in spatial.items():
	# 	w.writerow([key[0], key[1], val[0], val[1]]) 

	# w = csv.writer(open("umi_y.csv", "w"))
	# for key, val in yy.items():
	# 	w.writerow([key, val[0], val[1]]) 

	w = csv.writer(open("umi_x.csv", "w"))
	for key, val in xx.items():
		w.writerow([key, val[0], val[1]]) 

	# w = csv.writer(open("umi_pcr.csv", "w"))
	# for key, val in pcr.items():
	# 	w.writerow([key, val[0], val[1]])