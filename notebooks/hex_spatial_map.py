#import necessary packages
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import seaborn as sns
import numpy as np

import scanpy.api as sc







##define the hexagon plotting function to be called by other functions in file

def hexagon(x, y, c='k'):
	from matplotlib.patches import Polygon

	points = [[x + (2/3), y],
			  [x + (1/3), y + 0.5],
   			  [x - (1/3), y + 0.5],
              [x - (2/3), y],
   			  [x - (1/3), y - 0.5],
   			  [x + (1/3), y - 0.5]]
    
	polygon = Polygon(points, edgecolor='w', fc=c)
	return polygon

#in order to plot the hexagon patch: 
#plt.gca().add_patch(polygon)
#plt.axis('scaled')






#function for merging the matrix, annotation and gene annotation files and group cell data by spatial barcode.
#returns a count matrix with cell#, count, well barcode and gene name columns

def organize_matrix(matrix, cell_annotation, gene_annotation):
	import pandas as pd	

	#get and clean the cell annotation file
	annot=pd.read_csv(cell_annotation, header=None, names=['barcode','cell'])
	annot['barcode']=annot['barcode'].str[-16:] #seperation of barcode specific for current scheme


	#get the matrix file
	count_matrix=pd.read_csv(matrix, sep=' ', header=None, index_col=None, skiprows=3, names=['gene','cell','count'])

	#get the gene annotation file
	genes=pd.read_csv(gene_annotation, header=None, usecols=[3,4], names=['name','gene'])

	#merge the count matrix with the annotation file, then merge with the gene annotation file
	count_matrix_annotated=count_matrix.merge(annot, on='cell',how='left')

	count_matrix_annotated=count_matrix_annotated.merge(genes, on='gene',how='left')

	count_matrix_annotated.drop('gene',axis=1)

	return count_matrix_annotated









def gene_spatial_map(matrix, cell_annotation, gene_annotation, plate_map, map, cutoff):
	import pandas as pd

	#format the data and annotations
	count_matrix_annotated = organize_matrix(matrix, cell_annotation, gene_annotation)

	#get the barocode map
	spatial_map=pd.read_csv(plate_map, header=None, names=['barcode', 'X', 'Y'])


	#define the mapping type

	if map == 'UMIS':
		#group spatial barcodes and sum UMI counts per group
		count_matrix_annotated=count_matrix_annotated.groupby('barcode')['count'].sum().reset_index()


		map_array=spatial_map.merge(count_matrix_annotated, on='barcode', how='left')


		map_array=map_array.fillna('0')
		map_array['count']=map_array['count'].astype(int)



		#log transform
		map_array['count']=np.log(map_array['count']+1)

		#match colors to counts
		max_val=max(map_array['count'])

		norm = mpl.colors.Normalize(vmin=0,vmax=max_val)

		count_hex=[]

		#set threshold 6.124 = 500UMIs
		#              7.60 = 2000UMIs

		for i in map_array['count']:
			if i > cutoff:
 		   		count_hex.append(mpl.colors.rgb2hex(cm.hot(norm(i))[:3]))
			else:
				count_hex.append('#000000')

		map_array['color']=count_hex

		return map_array




	#elif: map = 'louvain':
		
	else:
		gene=map 
		
		gene_counts=pd.DataFrame(count_matrix_annotated.groupby(['barcode','name'])['count'].sum())


		specific_gene_counts=gene_counts.loc[pd.IndexSlice[:, gene], :]

		#drop the gene name index
		specific_gene_counts=specific_gene_counts.reset_index(level=1, drop=True)


		map_array=spatial_map.merge(specific_gene_counts, on='barcode', how='left')

		map_array=map_array.fillna('0')
		map_array['count']=map_array['count'].astype(int)



		#log transform
		map_array['count']=np.log(map_array['count']+1)

		#match colors to counts
		max_val=max(map_array['count'])

		norm = mpl.colors.Normalize(vmin=0,vmax=max_val)

		count_hex=[]

		#set threshold 6.124 = 500UMIs
		#              7.60 = 2000UMIs

		for i in map_array['count']:
			if i > cutoff:
 		   		count_hex.append(mpl.colors.rgb2hex(cm.viridis(norm(i))[:3]))
			else:
				count_hex.append('#000000')

		map_array['color']=count_hex

		return map_array













def plot_spatial_array(matrix, cell_annotation, gene_annotation, plate_map, map, cutoff):


	map_array=gene_spatial_map(matrix, cell_annotation, gene_annotation, plate_map, map, cutoff)

	gene=map

	plt.rcParams['figure.figsize'] = [15, 7.5]

	for index, row in map_array.iterrows():
		if row['X'] % 2 != 0:
			row['Y']=row['Y']-0.5
    
		polygon = hexagon(row['X'],row['Y'],row['color'])
		plt.gca().add_patch(polygon)


	plt.ylim(18,0) #flip y axis
	plt.axis('scaled')
	plt.title('Log Transformed, w/ cutoff') 
	plt.axis('off')

	plt.Figure()

#in order to plot output, run plt.show()































# #add section for cell number count per well
# #add colorbar
# def spatial_map_bulk(matrix, cutoff, cell_annotation, plate_map):
# #takes the output matrix file or scanpy anndata, cell annotation file a UMI cutoff and a plate map 


# 	####annotate the count matrix with spatial barcodes####

# 	#get and clean the cell annotation file
# 	cell_annotation_columns=['barcode','cell']
# 	annot=pd.read_csv(cell_annotation, header=None, names=cell_annotation_columns)
# 	annot['barcode']=annot['barcode'].str[-16:] #seperation of barcode specific for current scheme


# 	array_columns=['gene','cell','count']
# 	count_matrix=pd.read_csv(matrix, sep=' ', header=None, index_col=None, skiprows=3, names=array_columns)

# 	#merge the count matrix with the annotation file
# 	count_matrix_annotated=count_matrix.merge(annot, on='cell',how='left')

# 	#group spatial barcodes and sum UMI counts per group
# 	count_matrix_annotated=count_matrix_annotated.groupby('barcode')['count'].sum().reset_index()


# 	#get the barocode map
# 	map_columns=['barcode', 'X', 'Y']
# 	spatial_map=pd.read_csv(plate_map, header=None, names=map_columns)

# 	map_array=spatial_map.merge(count_matrix_annotated, on='barcode', how='left')


# 	map_array=map_array.fillna('0')
# 	map_array['count']=map_array['count'].astype(int)


# 	#log transform
# 	map_array['count']=np.log(map_array['count']+1)

# 	#match colors to counts
# 	max_val=max(map_array['count'])

# 	norm = mpl.colors.Normalize(vmin=0,vmax=max_val)

# 	count_hex=[]

# 	#set threshold 6.124 = 500UMIs
# 	#              7.60 = 2000UMIs

# 	for i in map_array['count']:
# 		if i > cutoff:
#  	   		count_hex.append(mpl.colors.rgb2hex(cm.hot(norm(i))[:3]))
# 		else:
# 			count_hex.append('#000000')

# 	map_array['color']=count_hex



# 	plt.rcParams['figure.figsize'] = [15, 7.5]

# #plotting chunck - common between functions

# 	for index, row in map_array.iterrows():
# 		if row['X'] % 2 != 0:
# 			row['Y']=row['Y']-0.5
    
# 		polygon = hexagon(row['X'],row['Y'],row['color'])
# 		plt.gca().add_patch(polygon)


# 	plt.ylim(18,0) #flip y axis
# 	plt.axis('scaled')
# 	plt.title('Log Transformed, w/ cutoff')
# 	plt.axis('off')

# 	plt.Figure()

#in order to plot output, run plt.show()







# def spatial_map(matrix, cutoff, cell_annotation, plate_map):
# #take in anndata after scanpy processing, 
# #genes, list of genes, louvain group etc (seperate function of louvain mapping?)
	























