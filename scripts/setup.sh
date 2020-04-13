#!/bin/bash --login
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50GB
#SBATCH --job-name setup
#SBATCH --output=%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Boleraceae_rnaseq/bin:$PATH"

#Set Variables
genome='ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/brassica_oleracea/dna/Brassica_oleracea.BOL.dna.toplevel.fa.gz'
annotations='ftp://ftp.ensemblgenomes.org/pub/plants/release-44/gff3/brassica_oleracea/Brassica_oleracea.BOL.44.gff3.gz'
samples=$(sed '1d' ../misc/samples.csv | cut -d ',' -f 1)

#Prepare reference
mkdir ref
cd ref
mkdir annotations STAR
echo "Downloading genome and annotations"

#Download genome
wget $genome -O TO1000.fa.gz
gunzip TO1000.fa.gz
samtools faidx TO1000.fa

#Download annotations
wget $annotations -O annotations/TO1000.gff.gz
gunzip annotations/TO1000.gff.gz

#Make gtf file
echo "Converting gff to gtf"
gffread annotations/TO1000.gff -T -o annotations/TO1000.gtf

#Create STAR index files
echo "Making index files"
STAR \
	--runThreadN 10 \
	--runMode genomeGenerate \
	--genomeDir ./STAR/ \
	--genomeFastaFiles TO1000.fa \
	--sjdbGTFfile annotations/TO1000.gtf \
	--sjdbOverhang 96

#Create sample directories, download, format data
cd ../
echo "Creating sample directories and downloading data"

for i in $samples
do
	echo "Sample $a"
	mkdir $i $i/fastq $i/alignment $i/job_reports
	cd $i/fastq
	#Download data from SRA
	sra_list=$(awk -v FS=',' -v a="$i" '{if ($1==a) print $2}' ../../../misc/samples.tsv)
	for x in $sra_list
	do
		echo "Downloading $x"
		prefetch --max-size 100000000 $x
		mv $i/*sra ./
		rmdir $i
		echo "Converting $i"
		fastq-dump --gzip --split-3 "$i".sra
	done
	cd ../../
done

echo "Finished setting up"

