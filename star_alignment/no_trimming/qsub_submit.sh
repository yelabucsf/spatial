#!/bin/bash                         

#$ -S /bin/bash                     
#$ -cwd                            
#$ -r y                            
#$ -j y                           
#$ -l mem_free=10G                 
#$ -l arch=linux-x64               
#$ -l netapp=2G,scratch=2G         
#$ -l h_rt=00:29:00 
#$ -t 1:44      


R2_file = sed "2{$SGEID};d" file_list.txt

/ye/netapp/jimmie.ye/tools/STAR/bin/Linux_x86_64/STAR --genomeDir /ye/yelabstore/10x.ref/refdata-cellranger-mm10-1.2.0/star --outSAMstrandField intronMotif --readFilesIn $R2_file --runThreadN 6
