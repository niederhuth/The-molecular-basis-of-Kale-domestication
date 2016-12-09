#!/bin/bash

echo "Starting"

#Prepare bowtie2 index
echo "Making bowtie2 index"
cd bowtie2
time /usr/local/samtools/1.2/samtools faidx Boleracea_v2.1.31.fa
time /usr/local/bowtie2/2.2.9/bowtie2-build Boleracea_v2.1.31.fa Boleracea_v2.1.31
cd ../

#Prepare gff file
echo "Format gff files for tophat and cufflinks"
cd misc/

cd ..

#Prepare transcriptome index
echo "Making transcriptome index"
cd transcriptome
export PATH=/usr/local/samtools/0.1.18/:${PATH}
export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
time /usr/local/tophat/2.0.14/bin/tophat --transcriptome-index Boleracea_v2.1.31
cd ../

echo "Done"
