#Define functions
#Function for making a table of results for two conditions from dds
makeResultsTable <- function(x,conditionA,conditionB,filter=FALSE){
    require(DESeq2)
    bml <- sapply(levels(dds$condition),function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl]))
    bml <- as.data.frame(bml)
    y <- results(x,contrast=c("condition",conditionA,conditionB),independentFiltering=filter)
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
  ggplot(x[x$fdr < cutoff,],aes(x=Significant/Annotated,y=reorder(Term,Significant/Annotated))) + 
    geom_point(aes(color=fdr,size=Significant)) + 
    theme_bw() +
    scale_color_continuous(low="red",high="blue") +
    xlab("Gene Ratio (DEGs/Annotated Genes)") + 
    ylab("") +
    labs(size="Gene Count",color="FDR") 
}
#Function for making combined dotplots of GOterms
GOdotplot2 <- function(x){
  ggplot(x) + 
    geom_point(aes(x=c("TvC"),y=Term,size=TvC_sig,color=TvC_FDR)) +
    geom_point(aes(x=c("TvK"),y=Term,size=TvK_sig,color=TvK_FDR)) +
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
    path <- c("../../figures_tables/goTerms/")
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
sampleTable <- read.csv("../../misc/sample_metadata.csv",header=T)
#Read in Counts Table
dds <- DESeqDataSetFromHTSeqCount(sampleTable,design= ~ condition)
#Prefilter
dds <- dds[rowSums(counts(dds)) > 1,]
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
pdf("../../figures_tables/sampleHeatMap.pdf",width=6,height=6,paper='special')
sampleHeatMap(rld)
dev.off()
#Make a PCA plot of samples
pdf("../../figures_tables/samplePCA.pdf",width=6,height=6,paper='special')
pcaPlot(rld)
dev.off()

#Make results tables for each pairwise comparison
resTvC <- makeResultsTable(dds,"TO1000","cabbage",filter=FALSE)
resTvK <- makeResultsTable(dds,"TO1000","kale",filter=FALSE)
resKvC <- makeResultsTable(dds,"kale","cabbage",filter=FALSE)
#Combine results tables
resfull <- as.data.frame(rbind(resTvC,resTvK,resKvC))
#Adjust p-values for all results
resfull$padj <- p.adjust(resfull$pval,method="BH")

#Extract raw counts and output a table
rawCounts <- data.frame(gene = row.names(counts(dds,normalized=FALSE)),
                        counts(dds,normalized=FALSE))
write.table(rawCounts, "../../figures_tables/raw_counts.tsv",sep="\t",
            quote=FALSE, row.names=TRUE)

#Extract and output a table of normalized counts
normalizedCounts <- counts(dds, normalized=TRUE)
all_genes <- data.frame(gene=row.names(normalizedCounts), normalizedCounts,
             as.data.frame(sapply(levels(dds$condition),
             function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl]))),
             resfull[resfull$sampleA=="TO1000" & resfull$sampleB=="cabbage",c(6,8)],
             resfull[resfull$sampleA=="TO1000" & resfull$sampleB=="kale",c(6,8)],
             resfull[resfull$sampleA=="kale" & resfull$sampleB=="cabbage",c(6,8)])
colnames(all_genes) <- c("gene","cabbage1", "cabbage2", "kale1", "kale2", "kale3",
                         "TO10001", "TO10002", "TO10003", "cabbage_mean",
                         "kale_mean", "TO1000_mean", "TvC_log2FC", "TvC_padj",
                         "TvK_log2FC", "TvK_padj", "KvC_log2FC", "KvC_padj" )
write.table(all_genes, "../../figures_tables/Gene_expression_table.tsv",sep="\t",
            quote=FALSE, row.names=FALSE)

#Extract significant DEGs
sig <- na.omit(resfull[resfull$padj <= 0.05 & resfull$log2FC >= 1 | resfull$padj <= 0.05 & resfull$log2FC <= -1,])
#Count number of DEGs for each comparison
table(sig$sampleA,sig$sampleB)
#For each pairwise comparison, extract the DEGs
TvCsig <- sig[sig$sampleA == "TO1000" & sig$sampleB == "cabbage",]
TvKsig <- sig[sig$sampleA == "TO1000" & sig$sampleB == "kale",]
KvCsig <- sig[sig$sampleA == "kale" & sig$sampleB == "cabbage",]

#Venn Diagram of DEGs
#Make a table of overlap between pairwise comparisons
d2 <- data.frame(id=unique(sig$id))
d2 <- data.frame(id=d2$id,TvC=ifelse(d2$id %in% sig[sig$sampleA == "TO1000" & sig$sampleB == "cabbage",]$id, 1, 0),
  TvK=ifelse(d2$id %in% sig[sig$sampleA == "TO1000" & sig$sampleB == "kale",]$id, 1, 0),
  KvC=ifelse(d2$id %in% sig[sig$sampleA == "kale" & sig$sampleB == "cabbage",]$id, 1, 0)
)
#Make Venn Diagram
library(VennDiagram)
pdf("../../figures_tables/comparisonVennDiagram.pdf",width=6,height=6,paper='special')
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

TvCgotermUP <- factor(as.integer(resTvC$id %in% TvCsig[TvCsig$log2FC > 1,]$id))
names(TvCgotermUP) <- resTvC$id
TvCgotermUP <- topGO(TvCgotermUP,goTerms,nodeSize=5,"TvC_up",writeData=TRUE)
TvCgotermDOWN <- factor(as.integer(resTvC$id %in% TvCsig[TvCsig$log2FC < -1,]$id))
names(TvCgotermDOWN) <- resTvC$id
TvCgotermDOWN <- topGO(TvCgotermDOWN,goTerms,nodeSize=5,"TvC_down",writeData=TRUE)

TvKgotermUP <- factor(as.integer(resTvK$id %in% TvKsig[TvKsig$log2FC > 1,]$id))
names(TvKgotermUP) <- resTvK$id
TvKgotermUP <- topGO(TvKgotermUP,goTerms,nodeSize=5,"TvK_up",writeData=TRUE)
TvKgotermDOWN <- factor(as.integer(resTvK$id %in% TvKsig[TvKsig$log2FC < -1,]$id))
names(TvKgotermDOWN) <- resTvK$id
TvKgotermDOWN <- topGO(TvKgotermDOWN,goTerms,nodeSize=5,"TvK_down",writeData=TRUE)

KvCgotermUP <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC > 1,]$id))
names(KvCgotermUP) <- resKvC$id
KvCgotermUP <- topGO(KvCgotermUP,goTerms,nodeSize=5,"KvC_up",writeData=TRUE)
KvCgotermDOWN <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC < -1,]$id))
names(KvCgotermDOWN) <- resKvC$id
KvCgotermDOWN <- topGO(KvCgotermDOWN,goTerms,nodeSize=5,"KvC_down",writeData=TRUE)

upGoterm <- merge(TvCgotermUP$BP,TvKgotermUP$BP,by.x="GO.ID",by.y="GO.ID")
upGoterm <- merge(upGoterm,KvCgotermUP$BP,by.x="GO.ID",by.y="GO.ID")
upSig <- upGoterm[upGoterm$fdr.x < 0.05 | upGoterm$fdr.y < 0.05 | upGoterm$fdr < 0.05, ]
upSig <- data.frame(Term = upSig$Term, TvC_FDR = upSig$fdr.x, TvC_sig = upSig$Significant.x, TvK_FDR = upSig$fdr.y, TvK_sig = upSig$Significant.y, KvC_FDR = upSig$fdr, KvC_sig = upSig$Significant)
upSig$TvC_FDR <- ifelse(upSig$TvC_FDR < 0.05, upSig$TvC_FDR, NA)
upSig$TvC_sig <- ifelse(upSig$TvC_FDR < 0.05, upSig$TvC_sig, NA)
upSig$TvK_FDR <- ifelse(upSig$TvK_FDR < 0.05, upSig$TvK_FDR, NA)
upSig$TvK_sig <- ifelse(upSig$TvK_FDR < 0.05, upSig$TvK_sig, NA)
upSig$KvC_FDR <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_FDR, NA)
upSig$KvC_sig <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_sig, NA)
ggsave("../../figures_tables/goTerms/Up_BP.pdf",plot=GOdotplot2(upSig))

upGoterm <- merge(TvCgotermUP$MF,TvKgotermUP$MF,by.x="GO.ID",by.y="GO.ID")
upGoterm <- merge(upGoterm,KvCgotermUP$MF,by.x="GO.ID",by.y="GO.ID")
upSig <- upGoterm[upGoterm$fdr.x < 0.05 | upGoterm$fdr.y < 0.05 | upGoterm$fdr < 0.05, ]
upSig <- data.frame(Term = upSig$Term, TvC_FDR = upSig$fdr.x, TvC_sig = upSig$Significant.x, TvK_FDR = upSig$fdr.y, TvK_sig = upSig$Significant.y, KvC_FDR = upSig$fdr, KvC_sig = upSig$Significant)
upSig$TvC_FDR <- ifelse(upSig$TvC_FDR < 0.05, upSig$TvC_FDR, NA)
upSig$TvC_sig <- ifelse(upSig$TvC_FDR < 0.05, upSig$TvC_sig, NA)
upSig$TvK_FDR <- ifelse(upSig$TvK_FDR < 0.05, upSig$TvK_FDR, NA)
upSig$TvK_sig <- ifelse(upSig$TvK_FDR < 0.05, upSig$TvK_sig, NA)
upSig$KvC_FDR <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_FDR, NA)
upSig$KvC_sig <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_sig, NA)
ggsave("../../figures_tables/goTerms/Up_MF.pdf",plot=GOdotplot2(upSig))

upGoterm <- merge(TvCgotermUP$CC,TvKgotermUP$CC,by.x="GO.ID",by.y="GO.ID")
upGoterm <- merge(upGoterm,KvCgotermUP$CC,by.x="GO.ID",by.y="GO.ID")
upSig <- upGoterm[upGoterm$fdr.x < 0.05 | upGoterm$fdr.y < 0.05 | upGoterm$fdr < 0.05, ]
upSig <- data.frame(Term = upSig$Term, TvC_FDR = upSig$fdr.x, TvC_sig = upSig$Significant.x, TvK_FDR = upSig$fdr.y, TvK_sig = upSig$Significant.y, KvC_FDR = upSig$fdr, KvC_sig = upSig$Significant)
upSig$TvC_FDR <- ifelse(upSig$TvC_FDR < 0.05, upSig$TvC_FDR, NA)
upSig$TvC_sig <- ifelse(upSig$TvC_FDR < 0.05, upSig$TvC_sig, NA)
upSig$TvK_FDR <- ifelse(upSig$TvK_FDR < 0.05, upSig$TvK_FDR, NA)
upSig$TvK_sig <- ifelse(upSig$TvK_FDR < 0.05, upSig$TvK_sig, NA)
upSig$KvC_FDR <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_FDR, NA)
upSig$KvC_sig <- ifelse(upSig$KvC_FDR < 0.05, upSig$KvC_sig, NA)
ggsave("../../figures_tables/goTerms/Up_CC.pdf",plot=GOdotplot2(upSig))

downGoterm <- merge(TvCgotermDOWN$BP,TvKgotermDOWN$BP,by.x="GO.ID",by.y="GO.ID")
downGoterm <- merge(downGoterm,KvCgotermDOWN$BP,by.x="GO.ID",by.y="GO.ID")
downSig <- downGoterm[downGoterm$fdr.x < 0.05 | downGoterm$fdr.y < 0.05 | downGoterm$fdr < 0.05, ]
downSig <- data.frame(Term = downSig$Term, TvC_FDR = downSig$fdr.x, TvC_sig = downSig$Significant.x, TvK_FDR = downSig$fdr.y, TvK_sig = downSig$Significant.y, KvC_FDR = downSig$fdr, KvC_sig = downSig$Significant)
downSig$TvC_FDR <- ifelse(downSig$TvC_FDR < 0.05, downSig$TvC_FDR, NA)
downSig$TvC_sig <- ifelse(downSig$TvC_FDR < 0.05, downSig$TvC_sig, NA)
downSig$TvK_FDR <- ifelse(downSig$TvK_FDR < 0.05, downSig$TvK_FDR, NA)
downSig$TvK_sig <- ifelse(downSig$TvK_FDR < 0.05, downSig$TvK_sig, NA)
downSig$KvC_FDR <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_FDR, NA)
downSig$KvC_sig <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_sig, NA)
ggsave("../../figures_tables/goTerms/Down_BP.pdf",plot=GOdotplot2(downSig))

downGoterm <- merge(TvCgotermDOWN$MF,TvKgotermDOWN$MF,by.x="GO.ID",by.y="GO.ID")
downGoterm <- merge(downGoterm,KvCgotermDOWN$MF,by.x="GO.ID",by.y="GO.ID")
downSig <- downGoterm[downGoterm$fdr.x < 0.05 | downGoterm$fdr.y < 0.05 | downGoterm$fdr < 0.05, ]
downSig <- data.frame(Term = downSig$Term, TvC_FDR = downSig$fdr.x, TvC_sig = downSig$Significant.x, TvK_FDR = downSig$fdr.y, TvK_sig = downSig$Significant.y, KvC_FDR = downSig$fdr, KvC_sig = downSig$Significant)
downSig$TvC_FDR <- ifelse(downSig$TvC_FDR < 0.05, downSig$TvC_FDR, NA)
downSig$TvC_sig <- ifelse(downSig$TvC_FDR < 0.05, downSig$TvC_sig, NA)
downSig$TvK_FDR <- ifelse(downSig$TvK_FDR < 0.05, downSig$TvK_FDR, NA)
downSig$TvK_sig <- ifelse(downSig$TvK_FDR < 0.05, downSig$TvK_sig, NA)
downSig$KvC_FDR <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_FDR, NA)
downSig$KvC_sig <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_sig, NA)
ggsave("../../figures_tables/goTerms/Down_MF.pdf",plot=GOdotplot2(downSig))

downGoterm <- merge(TvCgotermDOWN$CC,TvKgotermDOWN$CC,by.x="GO.ID",by.y="GO.ID")
downGoterm <- merge(downGoterm,KvCgotermDOWN$CC,by.x="GO.ID",by.y="GO.ID")
downSig <- downGoterm[downGoterm$fdr.x < 0.05 | downGoterm$fdr.y < 0.05 | downGoterm$fdr < 0.05, ]
downSig <- data.frame(Term = downSig$Term, TvC_FDR = downSig$fdr.x, TvC_sig = downSig$Significant.x, TvK_FDR = downSig$fdr.y, TvK_sig = downSig$Significant.y, KvC_FDR = downSig$fdr, KvC_sig = downSig$Significant)
downSig$TvC_FDR <- ifelse(downSig$TvC_FDR < 0.05, downSig$TvC_FDR, NA)
downSig$TvC_sig <- ifelse(downSig$TvC_FDR < 0.05, downSig$TvC_sig, NA)
downSig$TvK_FDR <- ifelse(downSig$TvK_FDR < 0.05, downSig$TvK_FDR, NA)
downSig$TvK_sig <- ifelse(downSig$TvK_FDR < 0.05, downSig$TvK_sig, NA)
downSig$KvC_FDR <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_FDR, NA)
downSig$KvC_sig <- ifelse(downSig$KvC_FDR < 0.05, downSig$KvC_sig, NA)
ggsave("../../figures_tables/goTerms/Down_CC.pdf",plot=GOdotplot2(downSig))


#KEGG analysis

#Genes of interest
path <- c("../../figures_tables/genes_of_interest/")
ifelse(!dir.exists(path),
dir.create(path), FALSE)
goi <- read.csv("../../misc/goi.csv")
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
