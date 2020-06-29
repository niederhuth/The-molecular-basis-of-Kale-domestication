TFs <- read.csv("../misc/Bo_TFs.csv",header=TRUE)

for(i in unique(TFs$TF.Family)){

}










tmp <- synRes[synRes$sampleA=="Kale" & synRes$sampleB=="TO1000",]
tmp$direction <- ifelse(tmp$log2FC > 0 & tmp$padj < 0.05,1,ifelse(tmp$log2FC < 0 & tmp$padj < 0.05,-1,0))
tmp2 <- tmp[tmp$direction != 0,]$At_gene
tmp2 <- as.vector(unique(tmp[tmp$direction != 0,]$At_gene))
tmp3 <- tmp[tmp$At_gene %in% tmp2,]

tmp4 <- tmp[,c(10,9,11)]
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$subgenome))
tmp5$sum <- rowSums(tmp5)
tmp6 <- as.vector(row.names(tmp5[tmp5$sum != tmp5$LF & tmp5$sum != tmp5$MF1 & tmp5$sum != tmp5$MF2,]))
tmp7 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$direction))
tmp7 <- tmp7[row.names(tmp5) %in% as.vector(unique(tmp$At_gene)),]
#tmp5 <- tmp5[row.names(tmp5) %in% tmp2,]
tmp7$sum <- rowSums(tmp7)

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




















