syn <- read.table("../misc/Bo_At_syntelogs.tsv",header=T,sep="\t")

nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="TO1000",])
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="TO1000",])
nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="TO1000",])/length(syn$Bo_gene %in% all_genes$gene)
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="TO1000",])/(nrow(all_genes)-length(syn$Bo_gene %in% all_genes$gene))
nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",])
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",])

nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="cabbage",])
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="cabbage",])
nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="cabbage",])/length(syn$Bo_gene %in% all_genes$gene)
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="cabbage",])/(nrow(all_genes)-length(syn$Bo_gene %in% all_genes$gene))
nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="cabbage",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",])
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="cabbage",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",])

nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])
nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])/length(syn$Bo_gene %in% all_genes$gene)
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])/(nrow(all_genes)-length(syn$Bo_gene %in% all_genes$gene))
nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])
nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])

nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF2',])

table(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id %in% syn[syn$subgenome=='LF',]$Bo_gene)[2]/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",])
table(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id %in% syn[syn$subgenome=='MF1',]$Bo_gene)[2]/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",])
table(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id %in% syn[syn$subgenome=='MF2',]$Bo_gene)[2]/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",])

nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='MF2',])

nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF2',])

nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF2',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="kale" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF2',])

nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='MF2',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="kale" & sig$sampleB=="cabbage",]$id,])/nrow(syn[syn$subgenome=='MF2',])

nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='LF' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='LF',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF1' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF1',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$log2FC>0 & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF2',])
nrow(syn[syn$subgenome=='MF2' & syn$Bo_gene %in% sig[sig$log2FC<0 & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]$id,])/nrow(syn[syn$subgenome=='MF2',])


