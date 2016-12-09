#!/bin/bash

echo "Starting"
mkdir raw_QC2 trimmed_QC2

#Unpack fastqs
cd fastq
bunzip2 *.fastq.bz2
cd ..

#QC analysis raw data
echo "First QC analysis"
export PATH=${PATH}:/usr/local/fastqc/latest/
time /usr/local/fastqc/latest/fastqc --noextract -t 10 -a ../../misc/evrogen_mint.tsv -o raw_QC2 fastq/*fastq 

#Trim data
echo "Trimming data"
time python2.7 /usr/local/cutadapt/1.9.dev1/bin/cutadapt --trim-n --max-n 0.1 -m 30 -q 20 -e 0.1 -O 3 -a file:../../misc/evogren_mint.fa  -o trimmed2.fastq fastq/*fastq

#QC analysis trimmed data
echo "Second QC analysis"
time /usr/local/fastqc/latest/fastqc --noextract -t 10 -a ../../misc/TruSeq.tsv -o trimmed_QC2 trimmed2.fastq 

#Align data
echo "Aligning data"
export PATH=/usr/local/samtools/0.1.18/:${PATH}
export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
time /usr/local/tophat/2.0.14/bin/tophat -F 0.1 --transcriptome-index ../ref/transcriptome/Gmax_v2 --no-coverage-search --b2-sensitive --no-mixed --read-realign-edit-dist 0 -i 30 -M -I 10000 -p 10 --library-type fr-firststrand -o tophat2 ../ref/bowtie2/Gmax_v2 trimmed2.fastq

#Repackage everything
#echo "Wrapping up"
#cd fastq
#for i in *fastq
#do
#bzip2 "$i"
#done
#cd ..

echo "Done"
