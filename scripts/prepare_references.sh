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
time /usr/local/cufflinks/2.2.1/bin/gffread Boleracea_v2.1.31.gff -T -o Boleracea_v2.1.31.gtf
time bash ../../../scripts/cuffdiff_gtf_attributes --input=Boleracea_v2.1.31.gtf --output=Boleracea_v2.1.31_fixed.gtf
cd ../

#Prepare transcriptome index
echo "Making transcriptome index"
cd transcriptome
export PATH=/usr/local/samtools/0.1.18/:${PATH}
export LD_LIBRARY_PATH=/usr/local/boost/1.54.0/gcc447/lib:/usr/local/gcc/4.7.1/lib:/usr/local/gcc/4.7.1/lib64:${LD_LIBRARY_PATH}
time /usr/local/tophat/2.0.14/bin/tophat -G ../misc/Boleracea_v2.1.31_fixed.gtf --transcriptome-index Boleracea_v2.1.31 ../bowtie2/Boleracea_v2.1.31
mv Boleracea_v2.1.31/* ../transcriptome/
rmdir Boleracea_v2.1.31/
rm -R tophat_out/
cd ../

echo "Done"
