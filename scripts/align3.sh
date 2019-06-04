#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB
#SBATCH --job-name align3
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/RNA/bin:$PATH"
#Define Variables
sample=$(pwd | sed s/.*\\///)
fastq='../fastq/trimmed.2.fastq.gz'

#map reads
cd alignment3
echo "Running STAR"
STAR \
 --runThreadN 20 \
 --runMode alignReads \
 --genomeDir ../../ref \
 --readFilesIn $fastq \
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
cut -f1,2 ReadsPerGene.out.tab | sed '1,4d' > counts.tsv

