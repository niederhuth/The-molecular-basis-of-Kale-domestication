---
Title: "The molecular basis of Kale domestication: Transcription profiling of developing leaves provides new insights into the evolution of a Brassica oleracea vegetative morphotype"
Authors: "Tatiana Arias, Chad E Niederhuth, Chris Pires, Paula McSteen"
Raw data: "GSE149483"
---
This repository is for scripts and processed data for the paper: 

Please cite this paper if you use any of the resources here.  

To reproduce the analysis, follow these steps:

**1) Clone this git repository**

```
git clone https://github.com/niederhuth/The-molecular-basis-of-kale-domestication-Comparative-transcriptomics
cd The-molecular-basis-of-kale-domestication-Comparative-transcriptomics
```

**NOTE:** These scripts were written for use on the MSU HPCC, which uses the SLURM workload manager. They also assume all software is installed in a conda environment. To run them on your computer or a different environment, you will need to change the header of each script. 

You will also need to modify this line in each script so that it points to your conda environment:

```
export PATH="$HOME/miniconda3/envs/Boleracea_rnaseq/bin:$PATH"
```

You _may_ also need to delete this line in each script:

```
cd $PBS_O_WORKDIR
```

**2) Create the conda environment**

```
conda env create -f scripts/Boleracea_rnaseq.yml
```

**3) Setup data directory and reference genome**

```
mkdir data
cd data
bash ../scripts/setup.sh
```
or submit as a job

```
mkdir data
cd data
sbatch ../scripts/setup.sh
```

**4) Trim reads**

```
cd <sample>
bash ../../scripts/trim.sh
```

or submit as a job

```
cd <sample>
sbatch ../../scripts/trim.sh
```

**5) Run the alignment**

This should be run in the data directory, above each of the samples. This is because the data are aligned in two passes. In there first pass, the data are aligned for each sample and novel junctions identified. These are grouped together for all samples and in the second pass each sample is realigned with the these novel junctions included.

```
cd ~/The-molecular-basis-of-kale-domestication-Comparative-transcriptomics/data
bash ../scripts/align.sh
```
or submit it as a job.

```
cd data
bash ../scripts/align.sh
```

This will produce the bam files and the <sample>\_counts.tsv that are the same as those found in the directory "figures_tables/raw_counts/"


**6) DEG analysis and figures** 

```
cd figures_tables
Rscript ../scripts/deg.R
```

To redo ReviGO, GO terms were pasted into: http://revigo.irb.hr/ and the Rscript downloaded to make the figures. These are all located in figures_tables/GO_terms/ReviGO_scripts/
You can run all these Rscripts using the ReviGO.R script

```
Rscript ../scripts/ReviGO.R
```

**Other Files**

* The **misc** directory contains a number of files used in this paper. This includes:

1) "adapters.fa" - adapter sequences
2) "Bo_annotations.csv" - descriptions of each _B. oleracea_ gene based on its Arabidopsis hit.
3) "Bo_At_syntelogs.csv" - classificaiton of _B. oleracea_ genes into its subgenome (LF, MF1, MF2) and its corresponding syntenic gene in Arabidopsis.
4) "Bo2ncbi.csv" - mapping of _B. oleracea_ genes to its NCBI accession for use in KEGG enrichment.
5) "deseq2_samples.csv" - input table describing experimental design for DESeq2.
6) "genes_of_interest.csv" - lists of particular genes we were interested in. Used for making heatmaps.
7) "samples.csv" - used in setting up directories (the setup.sh script) and downloading SRA files.
8) "topGO.txt" - GO term annotations for use with topGO in R.

* Other **scripts** files

1) "functions.R" - some R functions for frequently used steps. Used in the "deg.R" script.
2) "annotate_uniprot.sh" - downloads and runs blast against uniprot for GO term annotation. **Note** extracting the GO terms takes a lot more work, which I did manually. Just use the annotations I provide in the **misc** directory.
3) "annotate_ncbi.sh" - annotates this genome version with the NCBI accessions from the version in NCBI, needed for KEGG analysis. 


