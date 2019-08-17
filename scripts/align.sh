#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB
#SBATCH --job-name align
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Boleracea_rnaseq/bin:$PATH"

#Define variables
sample=$(pwd | sed s/.*\\///)
input="$sample.fastq.gz"
trimmed="trimmed.fastq.gz"
adaptors="file:../../../misc/adapters.fa"

#Cutadapt
cd fastq
echo "Running cutadapt"
cutadapt -j 20 --trim-n -m 30 -g $adaptors -o $trimmed $input

#Fastqc
mkdir fastqc
echo "Running fastqc"
fastqc -t 20 -o fastqc/ $input $trimmed

#map reads
cd ../alignment
echo "Running STAR"
STAR \
 --runThreadN 20 \
 --runMode alignReads \
 --genomeDir ../../ref \
 --readFilesIn $trimmed \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate \
 --outSAMstrandField intronMotif \
 --outFilterType BySJout \
 --outFilterMultimapNmax 20 \
 --alignSJoverhangMin 8 \
 --alignSJDBoverhangMin 1 \
 --alignIntronMin 20 \
 --alignIntronMax 1000000 \
 --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverReadLmax 0.04 \
 --quantMode GeneCounts

#make counts table
cut -f1,2 ReadsPerGene.out.tab | sed '1,4d' > "$sample"_counts.tsv

