#!/bin/bash

################################ pipeline parameters ###########################################

##############  these parameters will likely stay the same for each analysis  ##################

gtf_file="/ye/yelabstore2/spatial/genes.gtf" # file containing gene information to map to transcriptome 
script_path="/ye/yelabstore2/christa/py_scripts" # the script folder for called python scripts
core=6 # number of cores to be used 
cutoff=1  # number of unique reads cutoff for splitting single cells
mismatch=1 # number of mismatches allowed between R2 UMI observed barcode and those provided in $barcodes  
index="/ye/netapp/jimmie.ye/ref/refdata-cellranger-hg19_and_mm10-1.2.0/star" # STAR RNA-seq reference genomes 
python_path="/netapp/home/christacaggiano/netapp/home/christacaggiano/python/bin". # python path, should be python2
STAR="/ye/netapp/jimmie.ye/tools/STAR/bin/Linux_x86_64/STAR"
samtools_path="/ye/netapp/jimmie.ye/tools/samtools-1.3"

################################ experiment specific  parameters #################################
##########################  input/out parameters for each analysis   ##############################


sample_ID="/ye/yelabstore2/spatial/XYZ_20180122_456_2/sample_ID.txt"  # text file containing PCR barcode IDs 
barcodes="/ye/yelabstore2/spatial/sid/v1/XYZ_plate3_barcodes.csv" # the RT barcode list for splitting cells
fastq_folder="/ye/yelabstore2/spatial/XYZ_20180122_456_2/fastq" # the folder for fastq files
all_output_folder="/ye/yelabstore2/spatial/XYZ_20180122_456_2" # the output folder

##################################################################################################

