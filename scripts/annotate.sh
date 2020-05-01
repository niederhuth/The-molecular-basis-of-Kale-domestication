#!/bin/bash --login
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60GB
#SBATCH --job-name annotate
#SBATCH --output=job_reports/%x-%j.SLURMout

cd $PBS_O_WORKDIR
export PATH="$HOME/miniconda3/envs/phylo/bin:$PATH"

#Define variables
ftp="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz ftp://ftp.ensemblgenomes.org/pub/plants/release-43/fasta/brassica_oleracea/pep/Brassica_oleracea.BOL.pep.all.fa.gz"
uniprot_fa="uniprot_sprot.fasta"
brassica_fa="Brassica_oleracea.BOL.pep.all.fa"
uniprot_db="uniprot.dmnd"
output="blast.m8"
evalue=0.00001
max_target=1

#


#Download files and unzip
for i in $ftp
do
	wget $i
	gunzip $(echo $i | sed s/.*\\///)
done

#Make Diamond database
diamond makedb --db $uniprot_db --in $uniprot_fa

#BLAST Brassica oleracea against uniprot db
diamond blastp --max-target-seqs $max_target --evalue $evalue --db $uniprot_db --query $brassica_fa --out $output
