#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --job-name trim
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Boleracea_rnaseq/bin:$PATH"

#Define variables
sample=$(pwd | sed s/.*\\///)
input="$sample.fastq.gz"
output="trimmed.fastq.gz"
adaptors="file:../../../misc/adapters.fa"

#Cutadapt
cd fastq
echo "Running cutadapt"
cutadapt -j 10 --trim-n -m 30 -g $adaptors -o $output $input

#Fastqc
mkdir fastqc
echo "Running fastqc"
fastqc -t 10 -o fastqc/ $input $output

echo "Done"
