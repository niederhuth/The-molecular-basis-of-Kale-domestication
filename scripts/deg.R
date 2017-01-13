setwd("deg")

#Define functions
readCounts <- function() {
     l <- list.files()
     d <- read.table(l[1], header=F, sep="\t", row.names=1)
     for(i in l[2:length(l)]){
          d <- cbind(d, read.table(i, header=F, sep="\t", row.names=1))
     }
     colnames(d) <- gsub(pattern = "_counts.tsv", replacement = "", l)
     row.names(d) <- gsub(pattern = "gene:", replacement = "", row.names(d))
     return(d[1:(nrow(d)-5),])
}

sampleHeatMap <- function(x,y){
     require(gplots)
     require(RColorBrewer)
     dists = dist(t(exprs(x)))
     mat = as.matrix(dists)
     rownames(mat) = colnames(mat) = row.names(pData(y))
     hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
     heatmap.2(mat,trace="none",col = rev(hmcol),margin=c(13, 13),keysize=0.9,density.info="none",key.title="none",key.xlab="Sample distance")
}

sampleHeatMap2 <- function(x){
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

topGO <- function(genelist,goTerms,nodeSize,filename,writeData=FALSE){
    require(topGO)
    require(GO.db)
    path <- c("../../figures_tables/goTerms/")
    ifelse(!dir.exists(path),
    dir.create(path), FALSE)
    BP <- new("topGOdata",description="Biological Process",ontology="BP",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    MF <- new("topGOdata",description="Molecular Function",ontology="MF",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    CC <- new("topGOdata",description="Cellular Compartment",ontology="CC",
              allGenes=genelist,annot=annFUN.gene2GO,nodeSize=nodeSize,gene2GO=goTerms)
    FisherBP <- runTest(BP,algorithm="parentchild",statistic="fisher")
    FisherMF <- runTest(MF,algorithm="parentchild",statistic="fisher")
    FisherCC <- runTest(CC,algorithm="parentchild",statistic="fisher")
    BPgenTable <- GenTable(BP,Fisher=FisherBP,ranksOf="Fisher",topNodes=length(score(FisherBP)))
    MFgenTable <- GenTable(MF,Fisher=FisherMF,ranksOf="Fisher",topNodes=length(score(FisherMF)))
    CCgenTable <- GenTable(CC,Fisher=FisherCC,ranksOf="Fisher",topNodes=length(score(FisherCC)))
    BPgenTable$fdr <- p.adjust(BPgenTable$Fisher,method="BH")
    MFgenTable$fdr <- p.adjust(MFgenTable$Fisher,method="BH")
    CCgenTable$fdr <- p.adjust(CCgenTable$Fisher,method="BH")
    write.csv(BPgenTable,paste(path,filename,"_BP.csv",sep=""),row.names=FALSE,quote=FALSE)
    write.csv(MFgenTable,paste(path,filename,"_MF.csv",sep=""),row.names=FALSE,quote=FALSE)
    write.csv(CCgenTable,paste(path,filename,"_CC.csv",sep=""),row.names=FALSE,quote=FALSE)
    if(writeData){
      return(list(BP=BPgenTable,MF=MFgenTable,CC=CCgenTable))
    }
}


##
sampleTable <- read.csv("../../misc/sample_metadata.csv",header=T)
countTable <- readCounts()
condition <- factor(sampleTable$condition)

##
library(DESeq)
cds <- newCountDataSet(countTable,condition)
cds <- estimateSizeFactors(cds)

##Explore samples
cdsBlind <- estimateDispersions(cds,method="blind")
vsd <- varianceStabilizingTransformation(cdsBlind)
pdf("../../figures_tables/sampleHeatMap.pdf",width=6,height=6,paper='special')
sampleHeatMap(vsd,cdsBlind)
dev.off()

##Differential expression
cds <- estimateDispersions(cds,method="pooled",sharingMode="maximum",fitType="parametric")
TvC <- nbinomTest(cds,"TO1000", "cabbage")
TvC <- data.frame(TvC,sampleA=c("TO1000"),sampleB=c("cabbage"))
TvK <- nbinomTest(cds,"TO1000", "kale")
TvK <- data.frame(TvK,sampleA=c("TO1000"),sampleB=c("kale"))
KvC <- nbinomTest(cds,"kale", "cabbage")
KvC <- data.frame(KvC,sampleA=c("kale"),sampleB=c("cabbage"))
full <- as.data.frame(rbind(TvC,TvK,KvC))
full$padj <- p.adjust(full$pval,method="BH")

sig <- na.omit(full[full$padj <= 0.05 & full$log2FoldChange >= 2 | full$padj <= 0.05 & full$log2FoldChange <= -2,])
TCsig <- sig[sig$sampleA == "TO1000" & sig$sampleB == "cabbage",]
TKsig <- sig[sig$sampleA == "TO1000" & sig$sampleB == "kale",]
KCsig <- sig[sig$sampleA == "kale" & sig$sampleB == "cabbage",]
d <- data.frame(id=unique(sig$id))
d <- data.frame(id=d$id,TvC=ifelse(d$id %in% sig[sig$sampleA == "TO1000" & sig$sampleB == "cabbage",]$id, 1, 0),
  TvK=ifelse(d$id %in% sig[sig$sampleA == "TO1000" & sig$sampleB == "kale",]$id, 1, 0),
  KvC=ifelse(d$id %in% sig[sig$sampleA == "kale" & sig$sampleB == "cabbage",]$id, 1, 0)
)
library(VennDiagram)
pdf("../../figures_tables/comparisonVennDiagram.pdf",width=6,height=6,paper='special')
draw.triple.venn(area1=nrow(subset(d,TvC==1)),
  area2=nrow(subset(d,TvK==1)),
  area3=nrow(subset(d,KvC==1)),
  n12=nrow(subset(d,TvC==1 & TvK==1)),
  n23=nrow(subset(d,TvK==1 & KvC==1)),
  n13=nrow(subset(d,TvC==1 & KvC==1)),
  n123=nrow(subset(d,TvC==1 & TvK==1 & KvC==1)),
  category = c("TO1000 v Cabbage", "TO1000 v Kale", "Kale v Cabbage"),
  lty = "blank",
  fill = c("skyblue", "pink1", "mediumorchid")
)
dev.off()

##DESeq2
library(DESeq2)
dds <- DESeqDataSetFromHTSeqCount(sampleTable, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds,fitType="parametric")
dds <- nbinomWaldTest(dds)

##
rld <- rlog(dds, blind=FALSE)
rld <- rlog(dds, blind=TRUE)
pdf("../../figures_tables/sampleHeatMap2.pdf",width=6,height=6,paper='special')
sampleHeatMap2(rld)
dev.off()

pdf("../../figures_tables/samplePCA.pdf",width=6,height=6,paper='special')
pcaPlot(rld)
dev.off()

##Differential expression
resTvC <- makeResultsTable(dds,"TO1000","cabbage",filter=FALSE)
resTvK <- makeResultsTable(dds,"TO1000","kale",filter=FALSE)
resKvC <- makeResultsTable(dds,"kale","cabbage",filter=FALSE)
resfull <- as.data.frame(rbind(resTvC,resTvK,resKvC))
resfull$padj <- p.adjust(resfull$pval,method="BH")

sig2 <- na.omit(resfull[resfull$padj <= 0.05 & resfull$log2FoldChange >= 2 | resfull$padj <= 0.05 & resfull$log2FoldChange <= -2,])
table(sig2$sampleA,sig2$sampleB)
TCsig2 <- sig2[sig2$sampleA == "TO1000" & sig2$sampleB == "cabbage",]
TKsig2 <- sig2[sig2$sampleA == "TO1000" & sig2$sampleB == "kale",]
KCsig2 <- sig2[sig2$sampleA == "kale" & sig2$sampleB == "cabbage",]
d2 <- data.frame(id=unique(sig2$id))
d2 <- data.frame(id=d2$id,TvC=ifelse(d2$id %in% sig2[sig2$sampleA == "TO1000" & sig2$sampleB == "cabbage",]$id, 1, 0),
  TvK=ifelse(d2$id %in% sig2[sig2$sampleA == "TO1000" & sig2$sampleB == "kale",]$id, 1, 0),
  KvC=ifelse(d2$id %in% sig2[sig2$sampleA == "kale" & sig2$sampleB == "cabbage",]$id, 1, 0)
)

library(VennDiagram)
pdf("../../figures_tables/comparisonVennDiagram2.pdf",width=6,height=6,paper='special')
draw.triple.venn(area1=nrow(subset(d2,TvC==1)),
  area2=nrow(subset(d2,TvK==1)),
  area3=nrow(subset(d2,KvC==1)),
  n12=nrow(subset(d2,TvC==1 & TvK==1)),
  n23=nrow(subset(d2,TvK==1 & KvC==1)),
  n13=nrow(subset(d2,TvC==1 & KvC==1)),
  n123=nrow(subset(d2,TvC==1 & TvK==1 & KvC==1)),
  category = c("TO1000 v Cabbage", "TO1000 v Kale", "Kale v Cabbage"),
  lty = "blank",
  fill = c("skyblue", "pink1", "mediumorchid")
)
dev.off()

##GO term enrichment
library(topGO)
library(GO.db)
goTerms <- readMappings(file="../../misc/topGO.txt")

TvCgoterm <- factor(as.integer(resTvC$id %in% TCsig2$id))
names(TvCgoterm) <- resTvC$id
TvCgoterm <- topGO(TvCgoterm,goTerms,nodeSize=5,"TvC",writeData=TRUE)
TvKgoterm <- factor(as.integer(resTvK$id %in% TKsig2$id))
names(TvKgoterm) <- resTvK$id
TvKgoterm <- topGO(TvKgoterm,goTerms,nodeSize=5,"TvK",writeData=TRUE)
KvCgoterm <- factor(as.integer(resKvC$id %in% KCsig2$id))
names(KvCgoterm) <- resKvC$id
KvCgoterm <- topGO(KvCgoterm,goTerms,nodeSize=5,"KvC",writeData=TRUE)


ggplot(test,aes(x=c("TvC"),y=GO.ID)) + geom_point(aes(size=factor(Significant),color=factor(fdr)))
