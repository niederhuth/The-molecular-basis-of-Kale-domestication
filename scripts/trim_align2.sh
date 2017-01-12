#!/bin/bash

echo "Starting"
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)
mkdir trimmed_QC2

#Unpack fastqs
cd fastq
bunzip2 "$sample".fastq.bz2
cd ..

#Trim data
echo "Trimming data"
time python2.7 /usr/local/cutadapt/1.9.dev1/bin/cutadapt --quality-base=64 --trim-n --max-n 0.1 -m 30 -q 13,13 -e 0.1 -O 3 -a file:../../misc/evrogen_mint.fa  -o fastq/tmp.fastq fastq/"$sample".fastq
time /usr/local/fastx_toolkit/0.0.14/bin/fastx_trimmer -i fastq/tmp.fastq -o fastq/tmp2.fastq -f 2
perl /usr/local/prinseq/0.20.3/prinseq-lite.pl -verbose -fastq fastq/tmp2.fastq -min_len 30 -out_format 3 -out_good fastq/trimmed2 -out_bad fastq/tmp3
rm fastq/tmp.fastq fastq/tmp2.fastq fastq/tmp3.fastq

#QC analysis trimmed data
echo "Second QC analysis"
time /usr/local/fastqc/latest/fastqc --noextract -t 10 -a ../../misc/evrogen_mint.tsv -o trimmed_QC2 fastq/trimmed2.fastq

#Align data
echo "Aligning data"
export PATH=/usr/local/samtools/0.1.18/:${PATH}
export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
time /usr/local/tophat/2.0.14/bin/tophat -F 0.1 --phred64-quals --transcriptome-index ../ref/transcriptome/Boleracea_v2.1.31 --no-coverage-search --b2-sensitive --no-mixed --read-realign-edit-dist 0 -i 30 -M -I 10000 -p 10 --library-type fr-unstranded -o tophat2 ../ref/bowtie2/Boleracea_v2.1.31 fastq/trimmed2.fastq

#Count reads per feature
echo "Running htseq-count"
export PYTHONPATH=${PYTHONPATH}:/usr/local/htseq/0.6.1p1/lib/python
time /usr/local/htseq/0.6.1p1/bin/htseq-count -f bam -s no -t exon -i gene_id -m union tophat2/accepted_hits.bam ../ref/transcriptome/Boleracea_v2.1.31.gff > gene_counts2.tsv

#Repackage fastqs
echo "Wrapping up"
cd fastq
for i in *.fastq
bzip2 "$i"
done
cd ..

echo "Done"
