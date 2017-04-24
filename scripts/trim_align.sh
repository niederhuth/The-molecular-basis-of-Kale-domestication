#!/bin/bash

echo "Starting"
sample=$(pwd | sed s/.*data\\/// | sed s/\\/.*//)
mkdir fastq_QC

#Unpack fastqs
cd fastq
if [ -s trimmed2.fastq.gz ]
then
  gunzip trimmed2.fastq.gz
  cd ../
else
  for i in "$sample".fastq.gz
  do
    gunzip "$i"
  done

  #QC analysis raw data
  echo "First QC analysis"
  export PATH=${PATH}:/usr/local/fastqc/latest/
  time /usr/local/fastqc/0.11.4/fastqc --noextract -t 10 \
  -a ../../../misc/evrogen_mint.tsv -o ../raw_QC "$sample".fastq

  #Trim data
  echo "Trimming adapters"
  time python2.7 /usr/local/cutadapt/1.9.dev1/bin/cutadapt --quality-base=64 \
  --trim-n --max-n 0.1 -m 30 -q 20,20 -e 0.1 -O 3 \
  -a file:../../../misc/evrogen_mint.fa -o trimmed1.fastq "$sample".fastq

  #QC analysis adapter trimmed data
  echo "Second QC analysis"
  time /usr/local/fastqc/0.11.4/fastqc --noextract -t 10 \
  -a ../../../misc/evrogen_mint.tsv -o ../fastq_QC trimmed1.fastq

  #Second round of trimming
  echo "Second round of trimming"
  time perl /usr/local/prinseq/0.20.3/prinseq-lite.pl -verbose \
  -fastq tmp.fastq -trim_left 1 -trim_ns_right 1 -min_len 30 -out_format 3 \
  -out_good trimmed2 -out_bad tmp3
  rm trimmed1.fastq bad.fastq

  #QC analysis trimmed data
  echo "Third QC analysis"
  time /usr/local/fastqc/0.11.4/fastqc --noextract -t 10 \
  -a ../../../misc/evrogen_mint.tsv ../fastq_QC trimmed2.fastq
  cd ../
fi

#Align data
echo "Aligning data"
export PATH=/usr/local/samtools/0.1.18/:${PATH}
export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
time /usr/local/tophat/2.0.14/bin/tophat -F 0.1 --phred64-quals \
--transcriptome-index ../ref/transcriptome/Boleracea_v2.1.31 \
--no-coverage-search --b2-sensitive --no-mixed --read-realign-edit-dist 0 -i 30 \
-M -I 10000 -p 10 --library-type fr-unstranded \
-o tophat ../ref/bowtie2/Boleracea_v2.1.31 fastq/trimmed2.fastq

#Count reads per feature
echo "Counting reads mapping to genes"
export PYTHONPATH=${PYTHONPATH}:/usr/local/htseq/0.6.1p1/lib/python
time /usr/local/htseq/0.6.1p1/bin/htseq-count -f bam -s no -t exon -i gene_id \
-m union tophat/accepted_hits.bam ../ref/transcriptome/Boleracea_v2.1.31.gff \
> gene_counts.tsv

#Repackage everything
echo "Wrapping up"
cd fastq
for i in *fastq
do
  gzip "$i"
done
cd ..

echo "Done"
