# Spatial Seq Project @ the Ye Lab

### About

Pipeline for analyzing single cell RNA-seq data with spatial information. Adapted from pipeline @ the Trapnell lab

### Use

`qsub ./sciRNAseq.sh`

### Software Requirements

The following packages must be installed:
* [STAR](https://github.com/alexdobin/STAR) - Spliced Transcript Read Aligner  
* [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) - Read trimming software from Babraham Bioinformatics
* [SAMtools](http://samtools.sourceforge.net/) - suite of tools for interacting with sequencing data in SAM/BAM format
* [Python2](https://www.python.org/downloads/) - We know Python2 is kind of dumb, but this is just how it is right now

**Scripts** - /scripts/ - python path of scripts for this pipeline must be available for the pipeline

All required software paths should be put added to `parameters.sh`

Required parameters are:

```bash
# Pipeline parameters

$gtf_file="/gtf/path" # file path of gtf file
$script_path="folder/of/python/scripts" # location of python files to be used
$core=6 # number of cores to be used for experiment
$cutoff=1 # number of reads
$mismatch=1 # number of mismatches allowed between 'real'  and barcode and barcode observed in R2
$STAR="/path/to/star/executable/STAR" # path to STAR aligner
$trim_galore_path="/path/to/trimgalore/exe/TrimGalore" # path to trim galore
```

```bash
# Experiment parameters

$sample_ID="sample_id.txt" # file containing sample ID names
$barcodes="/folder/of/python/scripts" # location of python files to be used
$fastq_folder="/path/to/folder/with/fastq" # folder containing fastq.gz files
$all_output_folder="/output/folder" # output from file
```

To check the installation of these packages use `requirements.sh`
