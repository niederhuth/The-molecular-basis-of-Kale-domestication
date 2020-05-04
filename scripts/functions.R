#Define functions
#Function for making a table of results for two conditions from dds
makeResultsTable <- function(x,conditionA,conditionB,lfcThreshold=0,filter=FALSE){
    require(DESeq2)
    bml <- sapply(levels(dds$condition),function(lvl) rowMeans(counts(dds,
        normalized=TRUE)[,dds$condition == lvl]))
    bml <- as.data.frame(bml)
    y <- results(x,contrast=c("condition",conditionA,conditionB),
        lfcThreshold=lfcThreshold,independentFiltering=filter)
    y <- data.frame(id=gsub(pattern = "gene:",replacement = "",row.names(y)),
                    sampleA=c(conditionA),sampleB=c(conditionB),
                    baseMeanA=bml[,conditionA],baseMeanB=bml[,conditionB],
                    log2FC=y$log2FoldChange,pval=y$pvalue,padj=y$padj)
    row.names(y) <- c(1:nrow(y))
    return(y)
}


#Function for making a heat map of samples
sampleHeatMap <- function(x){
    require(pheatmap)
    require("RColorBrewer")
    sampleDists <- dist(t(assay(x)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(x)
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
}


#Function for making a PCA plot of samples
pcaPlot <- function(x){
    require(ggplot2)
    d <- plotPCA(x, intgroup=c("condition"), returnData=TRUE)
    pv <- round(100 * attr(d, "percentVar"))
    ggplot(d, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",pv[1],"% variance")) +
    ylab(paste0("PC2: ",pv[2],"% variance")) +
    coord_fixed() +
    theme(panel.background=element_blank(),
          axis.line=element_line(color="black"),
          axis.text=element_text(color="black"),
          axis.title=element_text(color="black",face="bold"))
}


#Function for making a heat map of genes
geneHeatMap <- function(dds,geneList){
    require(pheatmap)
    select <- select <- row.names(counts(dds,normalized=TRUE)) %in% genelist
    nt <- normTransform(dds) # defaults to log2(x+1)
    log2.norm.counts <- assay(nt)[select,]
    COL <- as.data.frame(colData(dds)[,c("condition")])
    pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=FALSE, annotation_col=df)
}


#Function for making dot plots of Enriched GO terms
GOdotplot <- function(x,fdr=0.05){
  require(ggplot2)
  x=head(x[x$fdr < fdr,],10)
  ggplot(x[x$fdr < fdr,],aes(x=Significant/Annotated,
    y=reorder(Term,Significant/Annotated))) + 
    geom_point(aes(color=fdr,size=Significant)) + 
    theme_bw() +
    scale_color_continuous(low="red",high="blue") +
    xlab("Gene Ratio (DEGs/Annotated Genes)") + 
    ylab("") +
    labs(size="Gene Count",color="FDR") +
    ggtitle("Top 10 Enriched GO Terms")
}
#Function for making combined dotplots of GOterms
GOdotplot2 <- function(x){
  ggplot(x) + 
    geom_point(aes(x=c("CvT"),y=Term,size=CvT_sig,color=CvT_FDR)) +
    geom_point(aes(x=c("KvT"),y=Term,size=KvT_sig,color=KvT_FDR)) +
    geom_point(aes(x=c("KvC"),y=Term,size=KvC_sig,color=KvC_FDR)) + 
    scale_color_continuous(low="red",high="blue",na.value="grey50") +
    theme_bw() +
    ylab("") +
    xlab("Sample Comparisons") +
    labs(size="Gene Count",color="FDR")
}


#Function for running topGO on a list of genes
topGO <- function(genelist,goTerms,nodeSize,fdr=0.05,filename,path="",writeData=FALSE){
    require(topGO)
    require(GO.db)
    ifelse(!dir.exists(path),dir.create(path), FALSE)
    BP <- new("topGOdata",description="Biological Process",ontology="BP",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    MF <- new("topGOdata",description="Molecular Function",ontology="MF",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    CC <- new("topGOdata",description="Cellular Compartment",ontology="CC",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    FisherBP <- runTest(BP,algorithm="parentchild",statistic="fisher")
    FisherMF <- runTest(MF,algorithm="parentchild",statistic="fisher")
    FisherCC <- runTest(CC,algorithm="parentchild",statistic="fisher")
    BPgenTable <- GenTable(BP,Fisher=FisherBP,ranksOf="Fisher",
        topNodes=length(score(FisherBP)))
    MFgenTable <- GenTable(MF,Fisher=FisherMF,ranksOf="Fisher",
        topNodes=length(score(FisherMF)))
    CCgenTable <- GenTable(CC,Fisher=FisherCC,ranksOf="Fisher",
        topNodes=length(score(FisherCC)))
    BPgenTable$fdr <- p.adjust(BPgenTable$Fisher,method="BH")
    MFgenTable$fdr <- p.adjust(MFgenTable$Fisher,method="BH")
    CCgenTable$fdr <- p.adjust(CCgenTable$Fisher,method="BH")
    write.csv(BPgenTable[BPgenTable$fdr <= fdr,],paste(path,filename,"_BP.csv",sep=""),
        row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_BP.pdf"),plot=GOdotplot(BPgenTable,fdr=fdr))
    write.csv(MFgenTable[MFgenTable$fdr <= fdr,],paste(path,filename,"_MF.csv",sep=""),
        row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_MF.pdf"),plot=GOdotplot(MFgenTable,fdr=fdr))
    write.csv(CCgenTable[CCgenTable$fdr <= fdr,],paste(path,filename,"_CC.csv",sep=""),
        row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_CC.pdf"),plot=GOdotplot(CCgenTable,fdr=fdr))
    if(writeData){
      return(list(BP=BPgenTable,MF=MFgenTable,CC=CCgenTable))
    }
}
