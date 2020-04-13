#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
#SBATCH --job-name trim
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Boleracea_rnaseq/bin:$PATH"

#Define variables
sample=$(pwd | sed s/.*\\///)
input=fastq/*.fastq.gz
trimmed="fastq/trimmed.fastq.gz"
adaptors="file:../../../misc/adapters.fa"
threads=20

#Cutadapt
echo "Running Cutadapt"
cutadapt \
	-j $threads \
	--trim-n \
	-m 30 \
	-g $adaptors \
	-o $trimmed $input

#Fastqc
mkdir fastq/fastqc
echo "Running fastqc"
fastqc \
	-t $threads \
	-o fastq/fastqc/ $input $trimmed

echo "Done"
