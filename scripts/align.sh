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
index="../../ref/STAR"
output1="rnaseq1"
output2="rnaseq2"
threads=20
fastq="../fastq/trimmed.fastq.gz"

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
		--outFilterMultimapNmax 10 \
		--alignSJoverhangMin 5 \
		--alignSJDBoverhangMin 3 \
		--alignIntronMin 20 \
		--alignIntronMax 0 \
		--outFilterScoreMinOverLread 0.33 \
		--outFilterMatchNminOverLread 0.33 \
		--outFilterMismatchNmax 10 \
		--outFilterMismatchNoverReadLmax 0.1
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
		--outFilterMultimapNmax 10 \
		--alignSJoverhangMin 5 \
		--alignSJDBoverhangMin 3 \
		--alignIntronMin 20 \
		--alignIntronMax 0 \
		--outFilterScoreMinOverLread 0.3 \
		--outFilterMatchNminOverLread 0.3 \
		--outFilterMismatchNmax 10 \
		--outFilterMismatchNoverReadLmax 0.1 \
		--quantMode GeneCounts
	cut -f1,4 ReadsPerGene.out.tab | sed '1,4d' > counts.tsv 
	cp counts.tsv ../../../figures_tables/raw_counts/"$i"_counts.tsv
	cd ../../
done

