#!/bin/bash --login
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20GB
#SBATCH --job-name align
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/Boleracea_rnaseq/bin:$PATH"

#Define variables
samples=$(sed '1d' ../misc/samples.csv | cut -d ',' -f 1 | tr '\n' ' ')
output1="rnaseq1"
output2="rnaseq2"
threads=20
fastq="../fastq/RNA/trimmed.fastq.gz"

#Star 1st pass
echo "Running star 1st pass"
for i in $samples
do
	echo "Running sample: $i"
	mkdir $i/$output1 $i/$output2
	cd $i/$output1
	STAR \
		--runThreadN $threads \
		--runMode alignReads \
		--genomeDir $index \
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
		--outFilterMismatchNoverReadLmax 0.04
	junctions="../../"$i"/"$output1"/SJ.out.tab $junctions"
	cd ../../
done
echo $junctions

#Star 2nd pass
echo "Running star 2nd pass"
for i in $samples
do
	echo "Running sample: $i"
	cd $i/$output2
	STAR \
		--runThreadN $threads \
		--runMode alignReads \
		--genomeDir $index \
		--readFilesIn "$fastq" \
		--sjdbFileChrStartEnd $junctions \
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
	cut -f1,4 ReadsPerGene.out.tab | sed '1,4d' > counts.tsv 
	cd ../../
done

