#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --job-name trim_reads
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/RNA/bin:$PATH"

#Sample
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)

#Input fastq file
p1="$sample.fastq.gz"

#Output fastq file
t1="trimmed.1.fastq.gz"

#Adapter sequences
#a1="AAGCAGTGGTATCAACGCAGAGTAC"
a1="TAGCAGTGGTATCAACGCAGAGTAC"

#Cutadapt
cd fastq
echo "Running cutadapt"
cutadapt -j 10 --trim-n -m 30 -a $a1 -o $t1 $p1

#Fastqc
mkdir fastqc
echo "Running fastqc"
fastqc -t 10 -o fastqc/ $t1 $p1

echo "Done"
