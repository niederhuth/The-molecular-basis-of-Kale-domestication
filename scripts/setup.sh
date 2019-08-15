#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60GB
#SBATCH --job-name setup
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Boleraceae_rnaseq/bin:$PATH"

#Create directories
echo "Setting things up"
mkdir data
cd data

#Download and prepare reference
mkdir ref
cd ref

echo "Downloading and genome and annotations"
#Download genome
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/brassica_oleracea/dna/Brassica_oleracea.BOL.dna.toplevel.fa.gz 
gunzip -c Brassica_oleracea.BOL.dna.toplevel.fa.gz > TO1000.fa
rm Brassica_oleracea.BOL.dna.toplevel.fa.gz
#Make fai
samtools faidx TO1000.fa
#Download annotations
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-44/gff3/brassica_oleracea/Brassica_oleracea.BOL.44.gff3.gz
gunzip -c Brassica_oleracea.BOL.44.gff3.gz > TO1000.gff
rm Brassica_oleracea.BOL.44.gff3.gz

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

#Create sample directories, download, format data
echo "Creating sample directories and downloading data"
cd ..

for a in cabbage1 cabbage2 kale1 kale2 kale3 TO10001 TO10002 TO10003
do
	echo "Sample $a"
	mkdir $a
	cd $a
	mkdir fastq alignment job_reports
	cd fastq
	#Download data from SRA
	python ../../../scripts/download_fastq.py $a
	for b in *sra
	do
		#Convert to fastq format
		fastq-dump --split-3
		cat *fastq | gzip > "$a".fastq.gz
		rm $b *fastq
	done
	cd ../../
done

echo "Finished setting up"

