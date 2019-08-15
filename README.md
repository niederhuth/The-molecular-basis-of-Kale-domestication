---
Title: "Analysis of _B. oleracea_ RNA-seq data."
RNA-seq Analysis: "Chad Niederhuth"
Collaborators: "Tatiana Arias, Chris Pires, Paula McSteen"
Raw data: "Link to SRA data will be posted at later date"
---
This repository is for scripts and processed data for the paper:

Please cite this paper if you use any of the resources here.  

All analyses performed on the Michigan State University High Performance Computing Cluster (HPCC)

To reproduce the analysis, follow these steps:
**NOTE:** These scripts were written for use on the MSU HPCC. To run them on your computer or a different environment, you will need to change the header of each script. You will also need to either modify or delete this line:
export PATH="$HOME/miniconda3/envs/Boleracea_rnaseq/bin:$PATH"
and delete this line
cd $PBS_O_WORKDIR

1) Clone this git repository

git clone https://github.com/niederhuth/The-molecular-basis-of-kale-domestication-Comparative-transcriptomics
cd The-molecular-basis-of-kale-domestication-Comparative-transcriptomics

2) Create the conda environment

conda env create -f Boleracea_rnaseq.yml

3) Run the setup script
bash scripts/setup.sh
or submit as a job

4) For each sample run the trim.sh script
cd <sample>
bash ../../scripts/trim.sh
or submit as a job

5) Run the alignment
cd <sample>
bash ../../scripts/align.sh
or submit as a job
This will produce the bam files and the <sample>_counts.tsv that are the same as those found in the directory "figures_tables/raw_counts/"

6) To rerun the DEG, functional analyses, and make figures
cd figures_tables
Rscript ../scripts/deg.R

**NOTE:** The annotate.sh script is solely for creating some of the GO terms and other annotation files found in the misc folder.



