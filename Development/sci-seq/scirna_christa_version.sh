#!/bin/bash
# this scRNA-seq pipeline accept a input folder, and then use the default parameter for the data processing and analysis, and generate a sparse gene count matrix for downstream analysis                                                                                                                                                                         
#                                                                                                                                                                                                                                                                                                                     
#$ -S /bin/bash                     
#$ -cwd                            
#$ -r y                            
#$ -j y                           
#$ -l mem_free=20G                 
#$ -l arch=linux-x64               
#$ -l netapp=2G,scratch=2G         
#$ -l h_rt=24:00:00

##sample_ID=$3 # the sample ID for each PCR samples after demultiplex
sample_ID=$1 # the sample ID for each PCR samples after demultiplex
##gtf_file=$4 # reference gtf files for gene counting
gtf_file=$2 # reference gtf files for gene counting
##core=$5 # core number for computation
core=$3 # core number for computation
##cutoff=$6  # the number of unique reads cutoff for splitting single cells
cutoff=$4  # the number of unique reads cutoff for splitting single cells
##index=$8 # STAR index for mapping
index="/ye/netapp/jimmie.ye/ref/star"
 
barcodes=/ye/yelabstore2/spatial/XYZ_20180122_cell_1_2/primary_barcodes_plate34.txt # the RT barcode list for splitting single cells
fastq_folder=/ye/yelabstore2/spatial/XYZ_20180122_cell_1_2/fastq # the folder for fastq files
all_output_folder=/ye/yelabstore2/spatial/XYZ_20180122_cell_1_2/ # the output folder
script_folder=/ye/yelabstore2/christa/py_scripts # the script folder for called python scripts

#define the mismatch rate (edit distance) of UMIs for removing duplicates:

mismatch=2

#define the bin of python
python_path="/netapp/home/christacaggiano/netapp/home/christacaggiano/python/bin"

#define the location of script:
script_path=$script_folder

# module load samtools/1.3
# module load bedtools/2.24.0

############ RT barcode and UMI attach
# this script take in a input folder, a sample ID, a output folder, a oligo-dT barcode file and
# call the python script to extract the UMI and RT barcode from read1 and attach them to the read names# of read2
input_folder=$fastq_folder
output_folder=$all_output_folder/UMI_attach
script=$script_folder/UMI_barcode_attach_gzipped.py
echo "changing the name of the fastq files..."

for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/$sample*R1*gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/$sample*R2*gz $input_folder/$sample.R2.fastq.gz; done

echo "Attaching barcode and UMI...."
mkdir -p $output_folder
$python_path/python2 $script $input_folder $sample_ID $output_folder $barcodes $core
echo "Barcode transformed and UMI attached."

# ################# Trim the read2
echo
echo "Start trimming the read2 file..."
echo $(date)
# module load python/2.7.3
module load cutadapt/1.8.3
module load trim_galore/0.4.1
mkdir $all_output_folder/trimmed_fastq
trimmed_fastq=$all_output_folder/trimmed_fastq
UMI_attached_R2=$all_output_folder/UMI_attach
for sample in $(cat $sample_ID); do echo trimming $sample; trim_galore $UMI_attached_R2/$sample*.gz -a AAAAAAAA --three_prime_clip_R1 1 -o $trimmed_fastq; done
# sem --semaphoretimeout 1800
echo "All trimmed file generated."
# module unload python/2.7.3

############align the reads with STAR, filter the reads, and remove duplicates based on UMI sequence and tagmentation site

#define the output folder
input_folder=$trimmed_fastq
STAR_output_folder=$all_output_folder/STAR_alignment
filtered_sam_folder=$all_output_folder/filtered_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam

#align read2 to the index file using STAR with default setting
echo "Start alignment using STAR..."
echo input folder: $input_folder
echo sample ID file: $sample_ID
echo index file: $index
echo output_folder: $STAR_output_folder

STAR="/ye/netapp/jimmie.ye/tools/STAR/bin/Linux_x86_64/STAR"
# make the output folder
mkdir -p $STAR_output_folder
# remove the index from the memory
# $STAR --genomeDir $index --genomeLoad Remove
# start the alignment
for sample in $(cat $sample_ID); do echo Aligning $sample;$STAR --runThreadN $core --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $input_folder/$sample*gz --outFileNamePrefix $STAR_output_folder/$sample --genomeLoad LoadAndKeep; done
# remove the index from the memory
# $STAR --genomeDir $index --genomeLoad Remove
echo "All alignment done."

# filter and sort the sam file
echo
echo "Start filter and sort the sam files..."
echo input folder: $STAR_output_folder
echo output folder: $filtered_sam_folder
mkdir -p $filtered_sam_folder

samtools_path="/ye/netapp/jimmie.ye/tools/samtools-1.3"
# module load samtools/1.3
for sample in $(cat $sample_ID); do echo Filtering $sample; $samtools_path/samtools view -bh -q 30 -F 4 $STAR_output_folder/$sample*.sam|$samtools_path/samtools sort -@ 10 -|$samtools_path/samtools view -h ->$filtered_sam_folder/$sample.sam; done

# Then for each filtered sam file, remove the duplicates based on UMI and barcode, chromatin number and position
echo
echo "Start removing duplicates..."
echo input folder: $filtered_sam_folder
echo output folder: $rmdup_sam_folder
mkdir -p $rmdup_sam_folder
# module unload python
for sample in $(cat $sample_ID); do echo remove duplicate $sample; $python_path/python2 $script_path/rm_dup_barcode_UMI.py $filtered_sam_folder/$sample.sam $rmdup_sam_folder/$sample.sam $mismatch; done 

# #mv the reported files to the report/duplicate_read/ folder
mkdir -p $input_folder/../report/duplicate_read
mv $rmdup_sam_folder/*.csv $input_folder/../report/duplicate_read/
mv $rmdup_sam_folder/*.png $input_folder/../report/duplicate_read/
echo "removing duplicates completed.."
echo
echo "Alignment and sam file preprocessing are done."  

################# split the sam file based on the barcode, and mv the result to the report folder
sam_folder=$all_output_folder/rmdup_sam
# bash $script_folder/samfile_split_multi_threads.sh $sam_folder $sample_ID $out_folder $barcodes $cutoff
sample_list=$sample_ID
output_folder=$all_output_folder/sam_splitted
barcode_file=$barcodes
cutoff=$cutoff

echo
echo "Start splitting the sam file..."
echo samfile folder: $sam_folder
echo sample list: $sample_list
echo ouput folder: $output_folder
echo barcode file: $barcode_file
echo cutoff value: $cutoff
mkdir -p $output_folder
# module unload python
for sample in $(cat $sample_list); do echo Now splitting $sample; $python_path/python2 $script_path/sam_split.py $sam_folder/$sample.sam $barcode_file $output_folder $cutoff; done


cat $output_folder/*sample_list.txt>$output_folder/All_samples.txt
cp $output_folder/All_samples.txt $output_folder/../barcode_samples.txt
# output the report the report/barcode_read_distribution folder
mkdir -p $output_folder/../report/barcode_read_distribution
mv $output_folder/*.txt $output_folder/../report/barcode_read_distribution/
mv $output_folder/*.png $output_folder/../report/barcode_read_distribution/
echo
echo "All sam file splitted."

################### calculate the reads number

fastq_folder=$fastq_folder
trimmed_folder=$trimmed_fastq
UMI_attach=$UMI_attached_R2
alignment=$STAR_output_folder
filtered_sam=$filtered_sam_folder
rm_dup_sam=$rmdup_sam_folder
#split_sam=$parental_folder/splited_sam
report_folder=$all_output_folder/report/read_num
echo
echo "Start calculating the reads number..."
#make the report folder
mkdir -p $report_folder
#calculate the read number and output the read number into the report folder
echo sample,total reads,after filtering barcode,after trimming,uniquely aligned reads,After remove duplicates>$report_folder/read_number.csv
for sample in $(cat $sample_ID); do echo calculating $sample; echo $sample,$(expr $(zcat $fastq_folder/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $UMI_attach/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $trimmed_folder/$sample*R2*.gz|wc -l) / 4),$($samtools_path/samtools view $filtered_sam/$sample.sam|wc -l),$(samtools view $rm_dup_sam/$sample.sam|wc -l)>>$report_folder/read_number.csv; done
echo "Read number calculation is done."

################## calculate the mouse and human and c.elegans reads fraction
input_folder=$all_output_folder/sam_splitted
sample_ID=$all_output_folder/barcode_samples.txt
output_folder=$all_output_folder/report/read_human_mouse
echo 
echo "Start calculating the mouse and human fraction..."
mkdir -p $output_folder
echo sample,human_reads,mouse_reads, cele_reads>$output_folder/human_mouse_fraction.txt
for sample in $(cat $sample_ID); do echo Processing $sample; echo $sample,$(samtools view $input_folder/$sample.sam|grep 'chr' -v|wc -l),$(samtools view $input_folder/$sample.sam|grep 'chr'|grep 'cele' -v|wc -l),$(samtools view $input_folder/$sample.sam|grep 'cele'|wc -l)>>$output_folder/human_mouse_fraction.txt; done
echo "Calculation done."

################# Generate the sparse gene count matrix
# count reads mapping to genes
output_folder=$all_output_folder/report/human_mouse_gene_count/
core_number=$core

script=$script_folder/sciRNAseq_count.py
echo "Start the gene count...."
$python_path/python2 $script $gtf_file $input_folder $sample_ID $core_number

echo "Make the output folder and transfer the files..."
mkdir -p $output_folder
# cat $input_folder/*.count > $output_folder/count.MM
# rm $input_folder/*.count
cat $input_folder/*.report > $output_folder/report.MM
rm $input_folder/*.report
mv $input_folder/*_annotate.txt $output_folder/
echo "All output files are transferred~"

echo "Analysis is done and gene count matrix is generated~"
