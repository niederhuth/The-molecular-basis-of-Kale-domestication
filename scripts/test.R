goi <- read.csv("../misc/goi.csv")

tmp <- synRes[synRes$sampleA=="kale" & synRes$sampleB=="TO1000",]
tmp$direction <- ifelse(tmp$log2FC > 0 & tmp$padj < 0.05,1,ifelse(tmp$log2FC < 0 & tmp$padj < 0.05,-1,0))
tmp2 <- tmp[tmp$direction != 0,]$At_gene
tmp2 <- as.vector(unique(tmp[tmp$direction != 0,]$At_gene))
tmp3 <- tmp[tmp$At_gene %in% tmp2,]

tmp4 <- tmp[,c(10,9,11)]
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$subgenome))
tmp5$sum <- rowSums(tmp5)
tmp6 <- as.vector(row.names(tmp5[tmp5$sum != tmp5$LF & tmp5$sum != tmp5$MF1 & tmp5$sum != tmp5$MF2,]))
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$direction))
tmp5 <- tmp5[row.names(tmp5) %in% as.vector(unique(tmp$At_gene)),]
#tmp5 <- tmp5[row.names(tmp5) %in% tmp2,]
tmp5$sum <- rowSums(tmp5)
tmp7 <- tmp5[tmp5$sum > 1 & row.names(tmp5) %in% tmp6,]
tmp6 <- row.names(tmp5[tmp5$`1` > 0 & tmp5$`-1` > 0,])
Bidirectional <- tmp[tmp$At_gene %in% tmp6,]
tmp6 <- row.names(tmp5[tmp5$sum == 1 & tmp5$`1` > 0,])
single_up <- tmp[tmp$At_gene %in% tmp6,]
tmp6 <- row.names(tmp5[tmp5$sum == 1 & tmp5$`-1` > 0,])
single_down <- tmp[tmp$At_gene %in% tmp6,]
tmp6 <- row.names(tmp5[tmp5$sum == 1,])
single_copy <- tmp[tmp$At_gene %in% tmp6,]

draw.triple.venn(area1=nrow(tmp7[tmp7$`-1`> 0,]),
                 area2=nrow(tmp7[tmp7$`1`> 0,]),
                 area3=nrow(tmp7[tmp7$`0`> 0,]),
                 n12=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0,]),
                 n23=nrow(tmp7[tmp7$`1`> 0 & tmp7$`0`> 0,]),
                 n13=nrow(tmp7[tmp7$`0`> 0 & tmp7$`-1`> 0,]),
                 n123=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0 & tmp7$`0`> 0,]),
                 category = c("Lower Expression", "Higher Expression", "No Difference"),
                 lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid")
)






nrow(tmp5[tmp5$sum > 1 & tmp5$`1` > 0 & tmp5$`-1` > 0,])


tmp[tmp$At_gene=="AT1G01050" & tmp$subgenome=="LF",]$direction


for(i in tmp$At_gene){
  tmp[tmp$At_gene==i,]
}





for(i in "Bo8g002920"){
  if(i %in% syn$Bo_gene){
    At_genes <- syn[syn$Bo_gene==i,]$At_gene
    Bo_genes <- syn[syn$At_gene==At_genes,]$Bo_gene
    #tmp <- synRes[synRes$id %in% Bo_genes,]
    tmp <- matrix(nrow=0,ncol=0)
    for(x in Bo_genes){
      if(x %in% resfull$id){
        tmp2 <- plotCounts(dds, x, intgroup = "condition", normalized = TRUE,transform = FALSE, returnData = TRUE)
        tmp2$gene <- x
        tmp <- rbind(tmp,tmp2)
      }
    }
    tmp <- merge(tmp,syn,by.x="gene",by.y="Bo_gene")
  }
  else{
    print("Not in syntenic genes")
  }
}
















pSyn$number=c(length(syn$Bo_gene %in% all_genes$gene),
          nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="TO1000",]),
          nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="cabbage",]),
          nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]),
          (nrow(all_genes)-length(syn$Bo_gene %in% all_genes$gene)),
          nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="TO1000",]),
          nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="cabbage",]),
          nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]))

test <- data.frame(row.names=c("Genome","KvT DEG","KvT Non-DEG","KvC DEG","KvC Non-DEG","CvT DEG","CvT Non-DEG"),
                   syntenic=c(24521,3519,21002,3469,21052,3514,21007),nonsyntenic=c(34704,4596,30108,3875,30829,3915,30789))


fisher.test(as.matrix(test[c(2:3),]))
fisher.test(as.matrix(test[c(4:5),]))
fisher.test(as.matrix(test[c(6:7),]))

fisher.test(as.matrix(test[c(2:3),c(2,1)]))
fisher.test(as.matrix(test[c(4:5),c(2,1)]))
fisher.test(as.matrix(test[c(6:7),c(2,1)]))

test2 <- data.frame(
  row.names=c("Genome","KvT DEG","KvT Non-DEG","KvC DEG","KvC Non-DEG","CvT DEG","CvT Non-DEG"),
  LF=c(10763,1463,9300,1472,9291,1502,9261),
  nonLF=c(48462,6652,41810,5872,42590,5526,42936),
  MF1=c(7420,1120,6300,1121,6299,1099,6321),
  nonMF1=c(51805,6995,44810,6223,45582,5929,45876),
  MF2=c(6338,936,5402,876,5462,913,5425),
  nonMF2=c(52887,7179,45708,6468,46419,6115,46772))

test3 <- data.frame(
  row.names=c("KvT","KvC","CvT"), 
  LF_pvalue=c(
    fisher.test(as.matrix(test2[c(2:3),c(1,2)]))$p.value,
    fisher.test(as.matrix(test2[c(4:5),c(1,2)]))$p.value,
    fisher.test(as.matrix(test2[c(6:7),c(1,2)]))$p.value), 
  LF_OR=c(
    fisher.test(as.matrix(test2[c(2:3),c(1,2)]))$estimate,
    fisher.test(as.matrix(test2[c(4:5),c(1,2)]))$estimate,
    fisher.test(as.matrix(test2[c(6:7),c(1,2)]))$estimate),
  MF1_pvalue=c(
    fisher.test(as.matrix(test2[c(2:3),c(3,4)]))$p.value,
    fisher.test(as.matrix(test2[c(4:5),c(3,4)]))$p.value,
    fisher.test(as.matrix(test2[c(6:7),c(3,4)]))$p.value), 
  MF1_OR=c(
    fisher.test(as.matrix(test2[c(2:3),c(3,4)]))$estimate,
    fisher.test(as.matrix(test2[c(4:5),c(3,4)]))$estimate,
    fisher.test(as.matrix(test2[c(6:7),c(3,4)]))$estimate),
  MF2_pvalue=c(
    fisher.test(as.matrix(test2[c(2:3),c(5,6)]))$p.value,
    fisher.test(as.matrix(test2[c(4:5),c(5,6)]))$p.value,
    fisher.test(as.matrix(test2[c(6:7),c(5,6)]))$p.value), 
  MF2_OR=c(
    fisher.test(as.matrix(test2[c(2:3),c(5,6)]))$estimate,
    fisher.test(as.matrix(test2[c(4:5),c(5,6)]))$estimate,
    fisher.test(as.matrix(test2[c(6:7),c(5,6)]))$estimate)
)

