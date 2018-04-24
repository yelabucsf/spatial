#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r y                            
#$ -j y                           
#$ -l mem_free=15G                 
#$ -l arch=linux-x64               
#$ -l netapp=2G,scratch=2G         
#$ -l h_rt=00:29:00 
#$ -t 1:44

# from a list of paths of read 2, get file corresponding to job #  
R2_file="`sed "${SGE_TASK_ID}q;d" file_list.txt`"

# echo $R2_file
# R2_file="/ye/yelabstore2/spatial/sid/v1/fastq_gz_all/A01.R2.fastq"

# make a directory for STAR output, and switch to that directory  
dir_name="A0$SGE_TASK_ID" 
mkdir -p "$dir_name"

trim_galore $R2_file --trim-n -a AAAAAAAA -clip_R1 9  -o $dir_name
cd "$dir_name" 

# run STAR with 6 threads. As in scRNA pipeline. 
/ye/netapp/jimmie.ye/tools/STAR/bin/Linux_x86_64/STAR --genomeDir /ye/yelabstore/10x.ref/refdata-cellranger-mm10-1.2.0/star --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outReadsUnmapped Fastx --outSAMstrandField intronMotif --readFilesCommand zcat --readFilesIn *.fq.gz --runThreadN 6
