#!/bin/bash

# SGE specific parameters:

#$ -S /bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=20G
#$ -l arch=linux-x64
#$ -l netapp=2G,scratch=2G
#$ -l h_rt=24:00:00

# single-cell RNA-seq pipeline adapted from C Trapnell's lab
# original pipeline here https://github.com/cole-trapnell-lab/single-cell-worm
# author <christa.caggiano@ucsf.edu> <sidraj@ucsf.edu> <jimmie.ye@ucsf.edu>


# source parameters that change with each run of pipeline
# parameters.sh must be in the same directory
source parameters.sh

# standardize naming of samples by naming them all R1 and R2
for sample in $(cat $sample_ID); do
	mv $input_folder/$sample*R1*gz $input_folder/$sample.R1.fastq.gz
	mv $input_folder/$sample*R2*gz $input_folder/$sample.R2.fastq.gz
done

############################## RT barcode and UMI attach #####################################

echo "Attaching barcode and UMI"
echo

output_folder=$all_output_folder/UMI_attach
script=$script_folder/UMI_barcode_attach_gzipped.py
mkdir -p $output_folder

# calls python script to attach UMI
if $python_path/python2 $script $input_folder $sample_ID $output_folder $barcodes $core then
	echo UMI assigned successfully
else
	echo ERROR: UMI assignment failed. Check fastq files.
	exit 1
fi
##################################### Trim R2  ###################################################

echo "Trimming R2"
echo

trimmed_fastq=$all_output_folder/trimmed_fastq
UMI_attached_R2=$all_output_folder/UMI_attach

mkdir $all_output_folder/trimmed_fastq

for sample in $(cat $sample_ID); do
	if $trim_galore_path $UMI_attached_R2/$sample*.gz -a AAAAAAAA --three_prime_clip_R1 1 -o $trimmed_fastq --suppress_warn then
		echo $sample trimmed successfully
	else
		echo ERROR: trimming of $sample failed. Please check UMI-attached files.
		exit 1
	fi
done

##################################### Align with STAR  ############################################

echo "Aligning with STAR"
echo

input_folder=$trimmed_fastq
STAR_output_folder=$all_output_folder/STAR_alignment
filtered_sam_folder=$all_output_folder/filtered_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam
STAR="/ye/netapp/jimmie.ye/tools/STAR/bin/Linux_x86_64/STAR"

mkdir -p $STAR_output_folder


for sample in $(cat $sample_ID); do

	if $STAR --genomeDir /ye/yelabstore2/christa/PHIX --outSAMstrandField intronMotif --readFilesCommand zcat
		--readFilesIn $input_folder/$sample*gz  --runThreadN 10 --outReadsUnmapped Fastx
		--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 then

		cp "Log.final.out" "Log.final.phix"

		$STAR --genomeDir $index --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
		--outSAMstrandField intronMotif  --readFilesIn Unmapped* --runThreadN 6

		cp "Log.final.out" "Log.final.mm10"

		echo $sample aligned successfully
	else
		echo ERROR: alignment failed. Please check input paramters.
	fi
done


##################################### Filter SAM files  #########################################

echo "Filtering sam files"
echo

samtools_path="/ye/netapp/jimmie.ye/tools/samtools-1.3"
mkdir -p $filtered_sam_folder

for sample in $(cat $sample_ID); do

	$samtools_path/samtools view -bh -q 30 -F 4 $STAR_output_folder/$sample*.sam|
	$samtools_path/samtools sort -@ 10 -|$samtools_path/samtools view -h ->$filtered_sam_folder/$sample.sam

done

##################################### Deduplication  #########################################

echo "Start removing duplicates..."
echo

mkdir -p $rmdup_sam_folder

for sample in $(cat $sample_ID); do

	if $python_path/python2 $script_path/rm_dup_barcode_UMI.py $filtered_sam_folder/$sample.sam $rmdup_sam_folder/$sample.sam $mismatch then
		echo $sample de-duplicated successfully.
	else
		echo ERROR: De-duplication of $sample failed. Please check sam file.
	fi
done

mkdir -p $input_folder/../report/duplicate_read
mv $rmdup_sam_folder/*.csv $input_folder/../report/duplicate_read/
mv $rmdup_sam_folder/*.png $input_folder/../report/duplicate_read/


############################## Split SAM files on barcode ###################################

echo "Splitting SAM file"
echo

sam_folder=$all_output_folder/rmdup_sam
sample_list=$sample_ID
output_folder=$all_output_folder/sam_splitted
barcode_file=$barcodes
cutoff=$cutoff
mkdir -p $output_folder

for sample in $(cat $sample_list); do

	if $python_path/python2 $script_path/sam_split.py $sam_folder/$sample.sam $barcode_file $output_folder $cutoff then
		echo $sample split correctly.
	else
		echo ERROR: check sam files. Splitting failed.
	fi
done


cat $output_folder/*sample_list.txt>$output_folder/All_samples.txt
cp $output_folder/All_samples.txt $output_folder/../barcode_samples.txt
mkdir -p $output_folder/../report/barcode_read_distribution
mv $output_folder/*.txt $output_folder/../report/barcode_read_distribution/
mv $output_folder/*.png $output_folder/../report/barcode_read_distribution/

############################## Calculate Reads ###################################

fastq_folder=$fastq_folder
trimmed_folder=$trimmed_fastq
UMI_attach=$UMI_attached_R2
alignment=$STAR_output_folder
filtered_sam=$filtered_sam_folder
rm_dup_sam=$rmdup_sam_folder
report_folder=$all_output_folder/report/read_num
mkdir -p $report_folder

echo "Calculating reads"

echo sample,total reads,after filtering barcode,after trimming,uniquely aligned reads,After remove duplicates>$report_folder/read_number.csv

for sample in $(cat $sample_ID); do

	echo $sample, $(expr $(zcat $fastq_folder/$sample*R2*.gz|wc -l) / 4),
	$(expr $(zcat $UMI_attach/$sample*R2*.gz|wc -l) / 4), $(expr $(zcat $trimmed_folder/$sample*R2*.gz|wc -l) / 4),
	$($samtools_path/samtools view $filtered_sam/$sample.sam|wc -l),$(samtools view $rm_dup_sam/$sample.sam|wc -l)>>$report_folder/read_number.csv;

done

############################## Calculate Human/Mouse fraction ###################################

echo "Calculating human and mouse fraction"
echo

input_folder=$all_output_folder/sam_splitted
sample_ID=$all_output_folder/barcode_samples.txt
output_folder=$all_output_folder/report/read_human_mouse
mkdir -p $output_folder

echo sample, human_reads, mouse_reads>$output_folder/human_mouse_fraction.txt
for sample in $(cat $sample_ID); do

	echo $sample, $(samtools view $input_folder/$sample.sam|grep 'chr' -v|wc -l),
	$(samtools view $input_folder/$sample.sam|grep 'chr'|grep 'cele' -v|wc -l),
	$(samtools view $input_folder/$sample.sam|grep 'cele'|wc -l)>>$output_folder/human_mouse_fraction.txt

done

############################## Generate gene count matrix ###################################


echo "Make gene count matrix"
echo

output_folder=$all_output_folder/report/human_mouse_gene_count/
core_number=$core
script=$script_folder/sciRNAseq_count.py


if $python_path/python2 $script $gtf_file $input_folder $sample_ID $core_number then
	echo genes counted successfully
else
	echo ERROR: gene counting failed
fi

mkdir -p $output_folder
cat $input_folder/*.count > $output_folder/count.MM
cat $input_folder/*.report > $output_folder/report.MM
rm $input_folder/*.report
mv $input_folder/*_annotate.txt $output_folder/

echo "Analysis is finished"
