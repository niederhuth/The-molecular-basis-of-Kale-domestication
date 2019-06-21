#Define functions
#Function for making a table of results for two conditions from dds
makeResultsTable <- function(x,conditionA,conditionB,lfcThreshold=0,filter=FALSE){
    require(DESeq2)
    bml <- sapply(levels(dds$condition),function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl]))
    bml <- as.data.frame(bml)
    y <- results(x,contrast=c("condition",conditionA,conditionB),lfcThreshold=lfcThreshold,independentFiltering=filter)
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
#Function for making dot plots of Enriched GO terms
GOdotplot <- function(x,cutoff=0.05){
  require(ggplot2)
  x=head(x[x$fdr < cutoff,],10)
  ggplot(x[x$fdr < cutoff,],aes(x=Significant/Annotated,y=reorder(Term,Significant/Annotated))) + 
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
topGO <- function(genelist,goTerms,nodeSize,filename,writeData=FALSE){
    require(topGO)
    require(GO.db)
    path <- c("goTerms/")
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
    BPgenTable <- GenTable(BP,Fisher=FisherBP,ranksOf="Fisher",topNodes=length(score(FisherBP)))
    MFgenTable <- GenTable(MF,Fisher=FisherMF,ranksOf="Fisher",topNodes=length(score(FisherMF)))
    CCgenTable <- GenTable(CC,Fisher=FisherCC,ranksOf="Fisher",topNodes=length(score(FisherCC)))
    BPgenTable$fdr <- p.adjust(BPgenTable$Fisher,method="BH")
    MFgenTable$fdr <- p.adjust(MFgenTable$Fisher,method="BH")
    CCgenTable$fdr <- p.adjust(CCgenTable$Fisher,method="BH")
    write.csv(BPgenTable,paste(path,filename,"_BP.csv",sep=""),row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_BP.pdf"),plot=GOdotplot(as.data.frame(BPgenTable)))
    write.csv(MFgenTable,paste(path,filename,"_MF.csv",sep=""),row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_MF.pdf"),plot=GOdotplot(MFgenTable))
    write.csv(CCgenTable,paste(path,filename,"_CC.csv",sep=""),row.names=FALSE,quote=FALSE)
    ggsave(paste(path,filename,"_CC.pdf"),plot=GOdotplot(CCgenTable))
    if(writeData){
      return(list(BP=BPgenTable,MF=MFgenTable,CC=CCgenTable))
    }
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

#DESeq2
#Run initial DESeq2 analysis
#Load Deseq2
library(DESeq2)
#Read in sample metadata
sampleTable <- read.csv("../misc/sample_metadata.csv",header=T)
#Read in Counts Table
dds <- DESeqDataSetFromHTSeqCount(sampleTable,design= ~ condition)
#Set reference level
dds$condition <- relevel(dds$condition, ref="TO1000")
#Estimate Size Factors
dds <- estimateSizeFactors(dds)
#Estimate Dispersions
dds <- estimateDispersions(dds,fitType="parametric")
#Run Wald Test
dds <- nbinomWaldTest(dds)

#Make diagnostic figures of samples
#Make a heat map of samples
rld <- rlog(dds, blind=TRUE)
pdf("sampleHeatMap.pdf",width=6,height=6,paper='special')
sampleHeatMap(rld)
dev.off()
#Make a PCA plot of samples
pdf("samplePCA.pdf",width=6,height=6,paper='special')
pcaPlot(rld)
dev.off()

#Make results tables for each pairwise comparison
resKvT <- makeResultsTable(dds,"kale","TO1000",lfcThreshold=0,filter=F)
resKvC <- makeResultsTable(dds,"kale","cabbage",lfcThreshold=0,filter=F)
resCvT <- makeResultsTable(dds,"cabbage","TO1000",lfcThreshold=0,filter=F)
#Combine results tables
resfull <- as.data.frame(rbind(resKvT,resKvC,resCvT))
#Adjust p-values for all results
resfull$padj <- p.adjust(resfull$pval,method="BH")

#Extract and output a table of normalized counts
normalizedCounts <- counts(dds, normalized=TRUE)
all_genes <- data.frame(gene=row.names(normalizedCounts), normalizedCounts[,c(3,4,5,1,2,6,7,8)],
             as.data.frame(sapply(levels(dds$condition),
             function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl]))),
             resfull[resfull$sampleA=="kale" & resfull$sampleB=="TO1000",c(6,8)],
             resfull[resfull$sampleA=="kale" & resfull$sampleB=="cabbage",c(6,8)],
             resfull[resfull$sampleA=="cabbage" & resfull$sampleB=="TO1000",c(6,8)])
colnames(all_genes) <- c("gene","kale1", "kale2", "kale3", "cabbage1", "cabbage2", 
                         "TO10001", "TO10002", "TO10003", "kale_mean",
                         "cabbage_mean", "TO1000_mean", "KvT_log2FC", "KvT_padj",
                         "KvC_log2FC", "KvC_padj", "CvT_log2FC", "CvT_padj" )
write.table(all_genes, "Gene_expression_table.tsv",sep="\t",
            quote=FALSE, row.names=FALSE)

#Extract significant DEGs
sig <- na.omit(resfull[resfull$padj <= 0.05 & resfull$log2FC >= 1 | resfull$padj <= 0.05 & resfull$log2FC <= -1,])
#Count number of DEGs for each comparison
table(sig$sampleA,sig$sampleB)
#For each pairwise comparison, extract the DEGs
CvTsig <- sig[sig$sampleA == "cabbage" & sig$sampleB == "TO1000",]
KvTsig <- sig[sig$sampleA == "kale" & sig$sampleB == "TO1000",]
KvCsig <- sig[sig$sampleA == "kale" & sig$sampleB == "cabbage",]

#Venn Diagram of DEGs
#Make a table of overlap between pairwise comparisons
d2 <- data.frame(id=unique(sig$id))
d2 <- data.frame(id=d2$id,KvT=ifelse(d2$id %in% sig[sig$sampleA == "kale" & sig$sampleB == "TO1000",]$id, 1, 0),
  KvC=ifelse(d2$id %in% sig[sig$sampleA == "kale" & sig$sampleB == "cabbage",]$id, 1, 0),
  CvT=ifelse(d2$id %in% sig[sig$sampleA == "cabbage" & sig$sampleB == "TO1000",]$id, 1, 0)
)
#Make Venn Diagram
library(VennDiagram)
pdf("comparisonVennDiagram.pdf",width=6,height=6,paper='special')
draw.triple.venn(area1=nrow(subset(d2,CvT==1)),
  area2=nrow(subset(d2,KvT==1)),
  area3=nrow(subset(d2,KvC==1)),
  n12=nrow(subset(d2,CvT==1 & KvT==1)),
  n23=nrow(subset(d2,KvT==1 & KvC==1)),
  n13=nrow(subset(d2,CvT==1 & KvC==1)),
  n123=nrow(subset(d2,CvT==1 & KvT==1 & KvC==1)),
  category = c("Cabbage v TO1000", "Kale v TO100", "Kale v Cabbage"),
  lty = "blank",
  fill = c("skyblue", "pink1", "mediumorchid")
)
dev.off()

##GO term enrichment
library(topGO)
library(GO.db)
goTerms <- readMappings(file="../misc/topGO.txt")

CvTgotermUP <- factor(as.integer(resCvT$id %in% CvTsig[CvTsig$log2FC > 1,]$id))
names(CvTgotermUP) <- resCvT$id
CvTgotermUP <- topGO(CvTgotermUP,goTerms,nodeSize=5,"CvT_up",writeData=TRUE)
CvTgotermDOWN <- factor(as.integer(resCvT$id %in% CvTsig[CvTsig$log2FC < -1,]$id))
names(CvTgotermDOWN) <- resCvT$id
CvTgotermDOWN <- topGO(CvTgotermDOWN,goTerms,nodeSize=5,"CvT_down",writeData=TRUE)

KvTgotermUP <- factor(as.integer(resKvT$id %in% KvTsig[KvTsig$log2FC > 1,]$id))
names(KvTgotermUP) <- resKvT$id
KvTgotermUP <- topGO(KvTgotermUP,goTerms,nodeSize=5,"KvT_up",writeData=TRUE)
KvTgotermDOWN <- factor(as.integer(resKvT$id %in% KvTsig[KvTsig$log2FC < -1,]$id))
names(KvTgotermDOWN) <- resKvT$id
KvTgotermDOWN <- topGO(KvTgotermDOWN,goTerms,nodeSize=5,"KvT_down",writeData=TRUE)

KvCgotermUP <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC > 1,]$id))
names(KvCgotermUP) <- resKvC$id
KvCgotermUP <- topGO(KvCgotermUP,goTerms,nodeSize=5,"KvC_up",writeData=TRUE)
KvCgotermDOWN <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC < -1,]$id))
names(KvCgotermDOWN) <- resKvC$id
KvCgotermDOWN <- topGO(KvCgotermDOWN,goTerms,nodeSize=5,"KvC_down",writeData=TRUE)

upGoterm <- merge(CvTgotermUP$BP,KvTgotermUP$BP,by.x="GO.ID",by.y="GO.ID")
upGoterm <- merge(upGoterm,KvCgotermUP$BP,by.x="GO.ID",by.y="GO.ID")
upSig <- upGoterm[upGoterm$fdr.x < 0.05 | upGoterm$fdr.y < 0.05 | upGoterm$fdr < 0.05, ]
upSig <- data.frame(Term = upSig$Term, CvT_FDR = upSig$fdr.x, CvT_sig = upSig$Significant.x, KvT_FDR = upSig$fdr.y, KvT_sig = upSig$Significant.y, KvC_FDR = upSig$fdr, KvC_sig = upSig$Significant)
upSig$CvT_FDR <- ifelse(upSig$CvT_FDR < 0.05, upSig$CvT_FDR, NA)
upSig$CvT_sig <- ifelse(upSig$CvT_FDR < 0.05, upSig$CvT_sig, NA)
upSig$KvT_FDR <- ifelse(upSig$KvT_FDR < 0.05, upSig$KvT_FDR, NA)
upSig$KvT_sig <- ifelse(upSig$KvT_FDR < 0.05, upSig$KvT_sig, NA)
upSig$KvC_FDR <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_FDR, NA)
upSig$KvC_sig <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_sig, NA)
ggsave("goTerms/Up_BP.pdf",plot=GOdotplot2(upSig))

upGoterm <- merge(CvTgotermUP$MF,KvTgotermUP$MF,by.x="GO.ID",by.y="GO.ID")
upGoterm <- merge(upGoterm,KvCgotermUP$MF,by.x="GO.ID",by.y="GO.ID")
upSig <- upGoterm[upGoterm$fdr.x < 0.05 | upGoterm$fdr.y < 0.05 | upGoterm$fdr < 0.05, ]
upSig <- data.frame(Term = upSig$Term, CvT_FDR = upSig$fdr.x, CvT_sig = upSig$Significant.x, KvT_FDR = upSig$fdr.y, KvT_sig = upSig$Significant.y, KvC_FDR = upSig$fdr, KvC_sig = upSig$Significant)
upSig$CvT_FDR <- ifelse(upSig$CvT_FDR < 0.05, upSig$CvT_FDR, NA)
upSig$CvT_sig <- ifelse(upSig$CvT_FDR < 0.05, upSig$CvT_sig, NA)
upSig$KvT_FDR <- ifelse(upSig$KvT_FDR < 0.05, upSig$KvT_FDR, NA)
upSig$KvT_sig <- ifelse(upSig$KvT_FDR < 0.05, upSig$KvT_sig, NA)
upSig$KvC_FDR <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_FDR, NA)
upSig$KvC_sig <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_sig, NA)
ggsave("goTerms/Up_MF.pdf",plot=GOdotplot2(upSig))

upGoterm <- merge(CvTgotermUP$CC,KvTgotermUP$CC,by.x="GO.ID",by.y="GO.ID")
upGoterm <- merge(upGoterm,KvCgotermUP$CC,by.x="GO.ID",by.y="GO.ID")
upSig <- upGoterm[upGoterm$fdr.x < 0.05 | upGoterm$fdr.y < 0.05 | upGoterm$fdr < 0.05, ]
upSig <- data.frame(Term = upSig$Term, CvT_FDR = upSig$fdr.x, CvT_sig = upSig$Significant.x, KvT_FDR = upSig$fdr.y, KvT_sig = upSig$Significant.y, KvC_FDR = upSig$fdr, KvC_sig = upSig$Significant)
upSig$CvT_FDR <- ifelse(upSig$CvT_FDR < 0.05, upSig$CvT_FDR, NA)
upSig$CvT_sig <- ifelse(upSig$CvT_FDR < 0.05, upSig$CvT_sig, NA)
upSig$KvT_FDR <- ifelse(upSig$KvT_FDR < 0.05, upSig$KvT_FDR, NA)
upSig$KvT_sig <- ifelse(upSig$KvT_FDR < 0.05, upSig$KvT_sig, NA)
upSig$KvC_FDR <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_FDR, NA)
upSig$KvC_sig <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_sig, NA)
ggsave("goTerms/Up_CC.pdf",plot=GOdotplot2(upSig))

downGoterm <- merge(CvTgotermDOWN$BP,KvTgotermDOWN$BP,by.x="GO.ID",by.y="GO.ID")
downGoterm <- merge(downGoterm,KvCgotermDOWN$BP,by.x="GO.ID",by.y="GO.ID")
downSig <- downGoterm[downGoterm$fdr.x < 0.05 | downGoterm$fdr.y < 0.05 | downGoterm$fdr < 0.05, ]
downSig <- data.frame(Term = downSig$Term, CvT_FDR = downSig$fdr.x, CvT_sig = downSig$Significant.x, KvT_FDR = downSig$fdr.y, KvT_sig = downSig$Significant.y, KvC_FDR = downSig$fdr, KvC_sig = downSig$Significant)
downSig$CvT_FDR <- ifelse(downSig$CvT_FDR < 0.05, downSig$CvT_FDR, NA)
downSig$CvT_sig <- ifelse(downSig$CvT_FDR < 0.05, downSig$CvT_sig, NA)
downSig$KvT_FDR <- ifelse(downSig$KvT_FDR < 0.05, downSig$KvT_FDR, NA)
downSig$KvT_sig <- ifelse(downSig$KvT_FDR < 0.05, downSig$KvT_sig, NA)
downSig$KvC_FDR <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_FDR, NA)
downSig$KvC_sig <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_sig, NA)
ggsave("goTerms/Down_BP.pdf",plot=GOdotplot2(downSig))

downGoterm <- merge(CvTgotermDOWN$MF,KvTgotermDOWN$MF,by.x="GO.ID",by.y="GO.ID")
downGoterm <- merge(downGoterm,KvCgotermDOWN$MF,by.x="GO.ID",by.y="GO.ID")
downSig <- downGoterm[downGoterm$fdr.x < 0.05 | downGoterm$fdr.y < 0.05 | downGoterm$fdr < 0.05, ]
downSig <- data.frame(Term = downSig$Term, CvT_FDR = downSig$fdr.x, CvT_sig = downSig$Significant.x, KvT_FDR = downSig$fdr.y, KvT_sig = downSig$Significant.y, KvC_FDR = downSig$fdr, KvC_sig = downSig$Significant)
downSig$CvT_FDR <- ifelse(downSig$CvT_FDR < 0.05, downSig$CvT_FDR, NA)
downSig$CvT_sig <- ifelse(downSig$CvT_FDR < 0.05, downSig$CvT_sig, NA)
downSig$KvT_FDR <- ifelse(downSig$KvT_FDR < 0.05, downSig$KvT_FDR, NA)
downSig$KvT_sig <- ifelse(downSig$KvT_FDR < 0.05, downSig$KvT_sig, NA)
downSig$KvC_FDR <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_FDR, NA)
downSig$KvC_sig <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_sig, NA)
ggsave("goTerms/Down_MF.pdf",plot=GOdotplot2(downSig))

downGoterm <- merge(CvTgotermDOWN$CC,KvTgotermDOWN$CC,by.x="GO.ID",by.y="GO.ID")
downGoterm <- merge(downGoterm,KvCgotermDOWN$CC,by.x="GO.ID",by.y="GO.ID")
downSig <- downGoterm[downGoterm$fdr.x < 0.05 | downGoterm$fdr.y < 0.05 | downGoterm$fdr < 0.05, ]
downSig <- data.frame(Term = downSig$Term, CvT_FDR = downSig$fdr.x, CvT_sig = downSig$Significant.x, KvT_FDR = downSig$fdr.y, KvT_sig = downSig$Significant.y, KvC_FDR = downSig$fdr, KvC_sig = downSig$Significant)
downSig$CvT_FDR <- ifelse(downSig$CvT_FDR < 0.05, downSig$CvT_FDR, NA)
downSig$CvT_sig <- ifelse(downSig$CvT_FDR < 0.05, downSig$CvT_sig, NA)
downSig$KvT_FDR <- ifelse(downSig$KvT_FDR < 0.05, downSig$KvT_FDR, NA)
downSig$KvT_sig <- ifelse(downSig$KvT_FDR < 0.05, downSig$KvT_sig, NA)
downSig$KvC_FDR <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_FDR, NA)
downSig$KvC_sig <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_sig, NA)
ggsave("goTerms/Down_CC.pdf",plot=GOdotplot2(downSig))

#KEGG analysis
library(clusterProfiler)
xx <- enrichMKEGG(KvTsig$id, organism='boe', minGSSize=1)


#Genes of interest
path <- c("genes_of_interest/")
ifelse(!dir.exists(path),
dir.create(path), FALSE)
goi <- read.csv("../misc/goi.csv")
for(gene in goi$gene){
    tryCatch({
        x <- plotCounts(dds, gene, intgroup = "condition",
                        normalized = TRUE,transform = FALSE, returnData = TRUE)
        p <- ggplot(x, aes(x = condition, y = count)) +
             geom_point(aes(color=condition)) +
             theme(panel.background=element_blank(),
                   axis.line=element_line(color="black"),
                   axis.text=element_text(color="black"),
                   axis.title=element_text(color="black",face="bold"),
                   legend.position="none") + xlab("Genotype")
    ggsave(paste(path,gene,"_counts.pdf",sep=""), p, width=5, height=4)}, error=function(e){})
}

#Syntenic genes
#
syn <- read.table("../misc/Bo-At-syntelogs.tsv",header=T,sep="\t")
synRes <- merge(resfull,syn,by.x="id",by.y="Bo_gene")
synSig <- merge(sig,syn,by.x="id",by.y="Bo_gene")

#Syntenic vs Non-Syntenic
pSyn=data.frame(
  row.names=c("Genome: Percent Syntenic",
              "KvT DEGs: Percent Syntenic",
              "KvC DEGs: Percent Syntenic",
              "CvT DEGs: Percent Syntenic",
              "Genome: Percent Non-Syntenic",
              "KvT DEGs: Percent Non-Syntenic",
              "KvC DEGs: Percent Non-Syntenic",
              "CvT DEGs: Percent Non-Syntenic"),
  percent=c(length(syn$Bo_gene %in% all_genes$gene)/nrow(all_genes),
            nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]),
            nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="kale" & sig$sampleB=="cabbage",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",]),
            nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]),
            (nrow(all_genes)-length(syn$Bo_gene %in% all_genes$gene))/nrow(all_genes),
            nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]),
            nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="kale" & sig$sampleB=="cabbage",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",]),
            nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])),
  order=c(1,2,3,4,5,6,7,8))
#Make Plot
ggsave("Syntenic_genes.pdf",
ggplot(pSyn,aes(x=reorder(row.names(pSyn),order),y=percent,fill=reorder(row.names(pSyn),order))) + 
  geom_bar(stat="identity") + 
  theme(panel.background=element_blank(),
        axis.line=element_line(color="black"),
        axis.text=element_text(color="black"),
        axis.title=element_text(color="black",face="bold"),
        axis.text.x=element_text(angle=315,hjust=0),
        legend.position="none") + 
  scale_y_continuous(expand=c(0,0),labels=percent) + 
  scale_fill_manual(values=c("tomato2",
                             "dodgerblue3",
                             "palegreen4",
                             "khaki3",
                             "tomato2",
                             "dodgerblue3",
                             "palegreen4",
                             "khaki3")) +
  ylab("Percentage of Genes") +
  xlab(""))

#KvT subgenome
KvTsub=data.frame(
  row.names=c("Genome: Percent LF",
              "DEGs: Percent LF",
              "Genome: Percent MF1",
              "DEGs: Percent MF1",
              "Genome: Percent MF2",
              "DEGs: Percent MF2"),
  percent=c(nrow(syn[syn$subgenome=="LF",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="LF",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF1",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF1",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF2",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF2",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="TO1000",])),
  order=c(1,2,3,4,5,6))
#Make Plot
ggsave("KvT_subgenome.pdf",
       ggplot(KvTsub,aes(x=reorder(row.names(KvTsub),order),y=percent,fill=reorder(row.names(KvTsub),order))) + 
         geom_bar(stat="identity") + 
         theme(panel.background=element_blank(),
               axis.line=element_line(color="black"),
               axis.text=element_text(color="black"),
               axis.title=element_text(color="black",face="bold"),
               axis.text.x=element_text(angle=315,hjust=0),
               legend.position="none") + 
         scale_y_continuous(expand=c(0,0),labels=percent) + 
         scale_fill_manual(values=c("tomato2",
                                    "tomato2",
                                    "dodgerblue3",
                                    "dodgerblue3",
                                    "black",
                                    "black")) +
         ylab("Percentage of Genes") +
         xlab(""))
#Log2FC distribution
tmp <- synRes[synRes$sampleA=="kale" & synRes$sampleB=="TO1000",]
tmp$direction <- ifelse(tmp$log2FC > 0 & tmp$padj < 0.05,1,ifelse(tmp$log2FC < 0 & tmp$padj < 0.05,-1,0))
tmp2 <- tmp[tmp$direction != 0,]$At_gene
tmp2 <- as.vector(unique(tmp[tmp$direction != 0,]$At_gene))
tmp3 <- tmp[tmp$At_gene %in% tmp2,]
ggsave("KvT_subgenome_log2FC.pdf",
  ggplot(tmp3[tmp3$padj<=0.05,],aes(x=log2FC,color=subgenome)) + 
    geom_density(size=1) +
    theme(panel.background=element_blank(),
          axis.line=element_line(color="black"),
          axis.text=element_text(color="black"),
          axis.title=element_text(color="black",face="bold"),
          axis.text.x=element_text(),
          legend.position="right") + 
    scale_y_continuous(expand=c(0,0)) +
    scale_color_manual(values=c("dodgerblue","tomato2","black")) +
    ylab("Density") +
    xlab("Log2 Fold Change"))
#Homoeolog direction
tmp4 <- tmp[,c(10,9,11)]
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$subgenome))
tmp5$sum <- rowSums(tmp5)
tmp6 <- as.vector(row.names(tmp5[tmp5$sum != tmp5$LF & tmp5$sum != tmp5$MF1 & tmp5$sum != tmp5$MF2,]))
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$direction))
tmp5 <- tmp5[row.names(tmp5) %in% as.vector(unique(tmp$At_gene)),]
tmp5$sum <- rowSums(tmp5)
tmp7 <- tmp5[tmp5$sum > 1 & row.names(tmp5) %in% tmp6,]
pdf("KvT_homoeologs_VennDiagram.pdf",width=6,height=6,paper='special')
draw.triple.venn(area1=nrow(tmp7[tmp7$`-1`> 0,]),
                 area2=nrow(tmp7[tmp7$`1`> 0,]),
                 area3=nrow(tmp7[tmp7$`0`> 0,]),
                 n12=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0,]),
                 n23=nrow(tmp7[tmp7$`1`> 0 & tmp7$`0`> 0,]),
                 n13=nrow(tmp7[tmp7$`0`> 0 & tmp7$`-1`> 0,]),
                 n123=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0 & tmp7$`0`> 0,]),
                 category = c("Lower Expression", "Higher Expression", "No Difference"),
                 lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

#KvC subgenome
KvCsub=data.frame(
  row.names=c("Genome: Percent LF",
              "DEGs: Percent LF",
              "Genome: Percent MF1",
              "DEGs: Percent MF1",
              "Genome: Percent MF2",
              "DEGs: Percent MF2"),
  percent=c(nrow(syn[syn$subgenome=="LF",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="kale" & synSig$sampleB=="cabbage" & synSig$subgenome=="LF",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",]),
            nrow(syn[syn$subgenome=="MF1",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="kale" & synSig$sampleB=="cabbage" & synSig$subgenome=="MF1",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",]),
            nrow(syn[syn$subgenome=="MF2",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="kale" & synSig$sampleB=="cabbage" & synSig$subgenome=="MF2",])/nrow(sig[sig$sampleA=="kale" & sig$sampleB=="cabbage",])),
  order=c(1,2,3,4,5,6))
#Make Plot
ggsave("KvC_subgenome.pdf",
       ggplot(KvCsub,aes(x=reorder(row.names(KvCsub),order),y=percent,fill=reorder(row.names(KvCsub),order))) + 
         geom_bar(stat="identity") + 
         theme(panel.background=element_blank(),
               axis.line=element_line(color="black"),
               axis.text=element_text(color="black"),
               axis.title=element_text(color="black",face="bold"),
               axis.text.x=element_text(angle=315,hjust=0),
               legend.position="none") + 
         scale_y_continuous(expand=c(0,0),labels=percent) + 
         scale_fill_manual(values=c("tomato2",
                                    "tomato2",
                                    "dodgerblue3",
                                    "dodgerblue3",
                                    "black",
                                    "black")) +
         ylab("Percentage of Genes") +
         xlab(""))
#Log2FC distribution
tmp <- synRes[synRes$sampleA=="kale" & synRes$sampleB=="cabbage",]
tmp$direction <- ifelse(tmp$log2FC > 0 & tmp$padj < 0.05,1,ifelse(tmp$log2FC < 0 & tmp$padj < 0.05,-1,0))
tmp2 <- tmp[tmp$direction != 0,]$At_gene
tmp2 <- as.vector(unique(tmp[tmp$direction != 0,]$At_gene))
tmp3 <- tmp[tmp$At_gene %in% tmp2,]
ggsave("KvC_subgenome_log2FC.pdf",
       ggplot(tmp3[tmp3$padj<=0.05,],aes(x=log2FC,color=subgenome)) + 
         geom_density(size=1) +
         theme(panel.background=element_blank(),
               axis.line=element_line(color="black"),
               axis.text=element_text(color="black"),
               axis.title=element_text(color="black",face="bold"),
               axis.text.x=element_text(),
               legend.position="right") + 
         scale_y_continuous(expand=c(0,0)) +
         scale_color_manual(values=c("dodgerblue","tomato2","black")) +
         ylab("Density") +
         xlab("Log2 Fold Change"))
#Homoeolog direction
tmp4 <- tmp[,c(10,9,11)]
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$subgenome))
tmp5$sum <- rowSums(tmp5)
tmp6 <- as.vector(row.names(tmp5[tmp5$sum != tmp5$LF & tmp5$sum != tmp5$MF1 & tmp5$sum != tmp5$MF2,]))
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$direction))
tmp5 <- tmp5[row.names(tmp5) %in% as.vector(unique(tmp$At_gene)),]
tmp5$sum <- rowSums(tmp5)
tmp7 <- tmp5[tmp5$sum > 1 & row.names(tmp5) %in% tmp6,]
pdf("KvC_homoeologs_VennDiagram.pdf",width=6,height=6,paper='special')
draw.triple.venn(area1=nrow(tmp7[tmp7$`-1`> 0,]),
                 area2=nrow(tmp7[tmp7$`1`> 0,]),
                 area3=nrow(tmp7[tmp7$`0`> 0,]),
                 n12=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0,]),
                 n23=nrow(tmp7[tmp7$`1`> 0 & tmp7$`0`> 0,]),
                 n13=nrow(tmp7[tmp7$`0`> 0 & tmp7$`-1`> 0,]),
                 n123=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0 & tmp7$`0`> 0,]),
                 category = c("Lower Expression", "Higher Expression", "No Difference"),
                 lty = "blank",
                 fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

#CvT subgenome
CvTsub=data.frame(
  row.names=c("Genome: Percent LF",
              "DEGs: Percent LF",
              "Genome: Percent MF1",
              "DEGs: Percent MF1",
              "Genome: Percent MF2",
              "DEGs: Percent MF2"),
  percent=c(nrow(syn[syn$subgenome=="LF",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="LF",])/nrow(sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF1",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF1",])/nrow(sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF2",])/nrow(all_genes),
            nrow(synSig[synSig$sampleA=="cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF2",])/nrow(sig[sig$sampleA=="cabbage" & sig$sampleB=="TO1000",])),
  order=c(1,2,3,4,5,6))
#Make Plot
ggsave("CvT_subgenome.pdf",
       ggplot(CvTsub,aes(x=reorder(row.names(CvTsub),order),y=percent,fill=reorder(row.names(CvTsub),order))) + 
         geom_bar(stat="identity") + 
         theme(panel.background=element_blank(),
               axis.line=element_line(color="black"),
               axis.text=element_text(color="black"),
               axis.title=element_text(color="black",face="bold"),
               axis.text.x=element_text(angle=315,hjust=0),
               legend.position="none") + 
         scale_y_continuous(expand=c(0,0),labels=percent) + 
         scale_fill_manual(values=c("tomato2",
                                    "tomato2",
                                    "dodgerblue3",
                                    "dodgerblue3",
                                    "black",
                                    "black")) +
         ylab("Percentage of Genes") +
         xlab(""))
#Log2FC distribution
tmp <- synRes[synRes$sampleA=="cabbage" & synRes$sampleB=="TO1000",]
tmp$direction <- ifelse(tmp$log2FC > 0 & tmp$padj < 0.05,1,ifelse(tmp$log2FC < 0 & tmp$padj < 0.05,-1,0))
tmp2 <- tmp[tmp$direction != 0,]$At_gene
tmp2 <- as.vector(unique(tmp[tmp$direction != 0,]$At_gene))
tmp3 <- tmp[tmp$At_gene %in% tmp2,]
ggsave("CvT_subgenome_log2FC.pdf",
       ggplot(tmp3[tmp3$padj<=0.05,],aes(x=log2FC,color=subgenome)) + 
         geom_density(size=1) +
         theme(panel.background=element_blank(),
               axis.line=element_line(color="black"),
               axis.text=element_text(color="black"),
               axis.title=element_text(color="black",face="bold"),
               axis.text.x=element_text(),
               legend.position="right") + 
         scale_y_continuous(expand=c(0,0)) +
         scale_color_manual(values=c("dodgerblue","tomato2","black")) +
         ylab("Density") +
         xlab("Log2 Fold Change"))
#Homoeolog direction
tmp4 <- tmp[,c(10,9,11)]
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$subgenome))
tmp5$sum <- rowSums(tmp5)
tmp6 <- as.vector(row.names(tmp5[tmp5$sum != tmp5$LF & tmp5$sum != tmp5$MF1 & tmp5$sum != tmp5$MF2,]))
tmp5 <- as.data.frame.matrix(table(tmp4$At_gene,tmp4$direction))
tmp5 <- tmp5[row.names(tmp5) %in% as.vector(unique(tmp$At_gene)),]
tmp5$sum <- rowSums(tmp5)
tmp7 <- tmp5[tmp5$sum > 1 & row.names(tmp5) %in% tmp6,]
pdf("CvT_homoeologs_VennDiagram.pdf",width=6,height=6,paper='special')
  draw.triple.venn(area1=nrow(tmp7[tmp7$`-1`> 0,]),
                  area2=nrow(tmp7[tmp7$`1`> 0,]),
                  area3=nrow(tmp7[tmp7$`0`> 0,]),
                  n12=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0,]),
                  n23=nrow(tmp7[tmp7$`1`> 0 & tmp7$`0`> 0,]),
                  n13=nrow(tmp7[tmp7$`0`> 0 & tmp7$`-1`> 0,]),
                  n123=nrow(tmp7[tmp7$`-1`> 0 & tmp7$`1`> 0 & tmp7$`0`> 0,]),
                  category = c("Lower Expression", "Higher Expression", "No Difference"),
                  lty = "blank",
                  fill = c("skyblue", "pink1", "mediumorchid"))
dev.off()

