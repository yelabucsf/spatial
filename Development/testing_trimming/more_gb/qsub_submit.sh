#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r y                            
#$ -j y                           
#$ -l mem_free=30G                 
#$ -l arch=linux-x64               
#$ -l netapp=2G,scratch=2G         
#$ -l h_rt=00:29:00 
#$ -t 4:7      


# from a list of paths of read 2, get file corresponding to job #  
R2_file="`sed "${SGE_TASK_ID}q;d" file_list.txt`"

echo $R2_file

# make a directory for STAR output, and switch to that directory  
dir_name="A0$SGE_TASK_ID" 
mkdir -p "$dir_name"
cd "$dir_name"

# run STAR with 6 threads. As in scRNA pipeline. 
/ye/netapp/jimmie.ye/tools/STAR/bin/Linux_x86_64/STAR --genomeDir /ye/yelabstore/10x.ref/refdata-cellranger-mm10-1.2.0/star --outSAMstrandField intronMotif --readFilesCommand zcat --readFilesIn $R2_file --runThreadN 6
