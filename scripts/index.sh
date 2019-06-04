#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --job-name index
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/RNA/bin:$PATH"

#Make gtf file
echo "Converting gff to gtf"
gffread TO1000.gff -T -o TO1000.gtf

#Create index files
echo "Making index files"
STAR \
 --runThreadN 10 \
 --runMode genomeGenerate \
 --genomeDir ./ \
 --genomeFastaFiles TO1000.fa \
 --sjdbGTFfile TO1000.gtf \
 --sjdbOverhang 96

