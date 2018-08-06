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

### Pipeline Output

``` bash
|── barcode_read_distribution
├── human_mouse_gene_count
│   ├── cell_annotate.txt
│   ├── count.MM
│   ├── gene_name_annotate.txt
│   ├── report.MM
│   └── report_annotate.txt
├── make_matrix.py
├── plate2_map.csv
├── read_human_mouse
│   └── human_mouse_fraction.txt
├── read_num
│   └── read_number.csv
└── spatial_output.csv
```

### summary files

* `make_matrix.py` makes `spatial_output.csv` summary file
* `spatial_output.csv` summary file that finds the number of transcripts for each well
    * Sample output:
    ``` bash
    # gene, x, y, PCR, number of counts
    8253,15,13,H10,3
  ```

#### human_mouse_gene_count
* `cell_annotate.txt` assigns an integer number to each cell (combination of spatial and PCR barcodes)  
    * Sample output:
    ``` bash
    # PCR.Spatial, cell_number
    A01.CTTACCGACAAGAGAA,1
    ```
* `gene_name_annotate.txt` assigns an integer number to a gene
    * Sample output:
    ``` bash
    # ensemble ID, function, exon/intron, name, gene #
    hg19_ENSG00000243485,"lincRNA",exon,hg19_MIR1302-10,1
    ```
* `report_annotate.txt` This file is a key for the categories of the type of gene annotations for each transcript on a scale of 1-8.
    * contents:
    ``` bash
    1, Perfect intersect exon match
    2, Nearest intersect exon match
    3, Perfect combine exon match
    4, Nearest combine exon match
    5, Perfect intersect gene match
    6, Nearest intersect gene match
    7, Perfect combine gene match
    8, Nearest combine gene match
    9, ambiguous match exons
    10, ambiguous match genes
    ```
* `report.MM` Assigns a category to each gene transcript, and the number of times that gene is observed
    * Sample output:
    ``` bash
    # gene #, report category, number of transcript observations
    1,8,664
    ```
* `count.MM` **most important file for downstream analysis** contains how many transcripts were observed for each cell/gene combination
    * Sample output:
    ``` bash
    # gene, cell number, number of observations
    8253,8,1
    ```

### Development

``` bash
├── barcode_maps
│   ├── XYZ_plate3_barcodes.csv
│   ├── calculate_map_intersection.py
│   ├── plate23_barcode_map.txt
│   ├── plate32_barcode_map.txt
│   ├── plate34_barcode_list.txt
│   ├── plate34_barcode_map.txt
│   ├── plate3_barcode_list.txt
│   ├── plate3_barcode_map.txt
│   ├── primary_barcodes_plate34.txt
│   └── primarybarcode34.txt
├── processing_scripts
│   ├── bin_human_mouse_by_spatial.py
│   ├── gene_matrix.py
│   ├── human-mouse.py
│   ├── jimmies-code.R
│   ├── percent_mapping.py
│   └── spatial_annotation_for_UMI_counts.py
├── sci-seq
│   ├── parameters.sh
│   ├── requirements.sh
│   ├── sciRNAseq.sh
│   └── scirna_christa_version.sh
```

* **barcode_maps/** contains all the barcode maps, if the experiment's map is not contained in its folder
* **processing_scripts/**
    * `bin_human_mouse_by_spatial.py`summarizes the number of human/mouse transcripts in a well **THIS SCRIPT IS A MESS**
    * `gene_matrix.py` makes a matrix of cell # by gene_ID
    * `human-mouse.py` sums the human/mouse read counts in various ways **ALSO MESS**
    * `jimmies-code.R` theoretically makes cool spatial plots, **pretty messy**
    * `percent_mapping.py` makes plot of percent mapping from STAR references
    * `spatial_annotation_for_UMI_counts.py` old version
* **sci-seq/** pipeline scripts 
