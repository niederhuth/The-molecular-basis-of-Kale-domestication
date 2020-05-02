evalue=0.00001
max_target=1
identity=70


wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/695/525/GCF_000695525.1_BOL/GCF_000695525.1_BOL_translated_cds.faa.gz

gunzip GCF_000695525.1_BOL_translated_cds.faa.gz 

sed s/\>.*GeneID\:/\>/ GCF_000695525.1_BOL_translated_cds.faa | \
	sed s/\].*// | fasta_formatter -w 0 > cds.faa

grep \> cds.faa | sort | uniq | while read line
do
	grep -m 1 -A 1 "$line" cds.faa >> cds2.faa
done

diamond makedb --db tcds.dmnd --in cds2.faa 
diamond makedb --db Bo.dmnd --in Brassica_oleracea.BOL.pep.all.fa

diamond blastp \
	--max-target-seqs $max_target \
	--evalue $evalue \
	--id $identity \
	--db Bo.dmnd \
	--query cds2.faa \
	--out tcds.m8

diamond blastp \
	--max-target-seqs $max_target \
	--evalue $evalue \
	--id $identity \
	--db tcds.dmnd \
	--query Brassica_oleracea.BOL.pep.all.fa \
	--out tcds2.m8

cut -f1 tcds2.m8 > tmp
cut -f2 tcds2.m8 > tmp2
paste tmp2 tmp > tmp3
fgrep -f tmp3 tcds.m8 > tmp4
cut -f1,2 tmp4 | sort | uniq > tmp5



