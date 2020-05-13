#Import user-defined functions
source("../scripts/functions.R")
library(Cairo)

#Create a directory for finalized paper figures & tables
path <- "Paper/"
if(!dir.exists(path)){
	dir.create(path)
}
#Make pdf table of mapping statistics
library(scales)
library(flextable)
mapStats <- read.csv("Mapping_Stats.csv",header=TRUE)
#Add a percentage mapped column
mapStats$Pecent.Uniquely.Mapped <- percent(
  mapStats$Uniquely.Mapped/mapStats$Trimmed.Reads,accuracy=0.1)
#Get rid of period in the column names
colnames(mapStats) <- gsub("\\."," ",colnames(mapStats))
#Create table with flextable
mapStatsFT <- flextable(mapStats)
mapStatsFT <- autofit(mapStatsFT)
mapStatsFT <- align(mapStatsFT,align="center",part="all")
mapStatsFT <- add_header_lines(mapStatsFT,
  values=("Supplementary Table 1: Mapping Statistics"))
mapStatsFT <- bold(mapStatsFT,part="header")
mapStatsFT <- bold(mapStatsFT,j="Sample")
mapStatsFT <- font(mapStatsFT,fontname="Arial")
mapStatsFT <- fontsize(mapStatsFT,size=12,part="all")
save_as_docx(mapStatsFT,path=paste(path,"Supplementary_Table_1.docx",sep=""))


#Run DESeq2 analysis
#Load Deseq2
library(DESeq2)
#Read in sample table
sampleTable <- read.csv("../misc/deseq2_samples.csv",header=T)
#Read in counts table
dds <- DESeqDataSetFromHTSeqCount(sampleTable,design= ~ condition)
#Prefilter data (optional)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#Set reference level
dds$condition <- relevel(dds$condition, ref="TO1000")
#Estimate size factors
dds <- estimateSizeFactors(dds)
#Estimate dispersions
dds <- estimateDispersions(dds,fitType="parametric")
#Run Wald test
dds <- nbinomWaldTest(dds)


#Make diagnostic figures of samples
path <- "Diagnostics/"
if(!dir.exists(path)){
  dir.create(path)
}
#Make a heat map of samples
rld <- rlog(dds, blind=TRUE)
cairo_pdf(paste(path,"Sample_Comparison_Heat_Map.pdf",sep=""),width=6,height=6,
  family="Arial")
	#sampleHeatMap is a user-defined function
	sampleHeatMap(rld)
dev.off()
#Make a PCA plot of samples
cairo_pdf(paste(path,"Sample_Comparison_PCA.pdf",sep=""),width=6,height=6,
  family="Arial")
	#pcaPlot is a user-defined function
	pcaPlot(rld) + guides(color=guide_legend(title="Genotype"))
dev.off()

#DEGs
#Make directory for DEG analyses
path <- "DEGs/"
if(!dir.exists(path)){
  dir.create(path)
}
#Make results tables for each pairwise comparison
#makeResultsTable is a user-defined function
resKvT <- makeResultsTable(dds,"Kale","TO1000",lfcThreshold=0,filter=F)
resKvC <- makeResultsTable(dds,"Kale","Cabbage",lfcThreshold=0,filter=F)
resCvT <- makeResultsTable(dds,"Cabbage","TO1000",lfcThreshold=0,filter=F)
#Combine results tables
resfull <- as.data.frame(rbind(resKvT,resKvC,resCvT))
#Adjust p-values for all results
#This is necessary because the adjusted p-value is only for the 
#individual pairwise comparisons
resfull$padj <- p.adjust(resfull$pval,method="BH")


#Extract and output a table of normalized counts
normalizedCounts <- counts(dds, normalized=TRUE)
geneEx <- data.frame(
				Gene=row.names(normalizedCounts),normalizedCounts[,c(3,4,5,1,2,6,7,8)],
				as.data.frame(sapply(levels(dds$condition),
					function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$condition == lvl]))),
				resfull[resfull$sampleA=="Kale" & resfull$sampleB=="TO1000",c(6,8)],
				resfull[resfull$sampleA=="Kale" & resfull$sampleB=="Cabbage",c(6,8)],
				resfull[resfull$sampleA=="Cabbage" & resfull$sampleB=="TO1000",c(6,8)]
				)
#Reset the column names
colnames(geneEx) <- c("Gene","Kale_1","Kale_2","Kale_3",
	"Cabbage_1","Cabbage_2","TO1000_1","TO1000_2","TO1000_3",
	"Kale_Mean","Cabbage_Mean","TO1000_Mean","KvT_log2FC",
	"KvT_padj","KvC_log2FC","KvC_padj","CvT_log2FC","CvT_padj")
#Output the table
write.table(geneEx,paste(path,"Gene_Expression_Table.csv",sep=""),sep=",",
  quote=FALSE,row.names=FALSE)


#Extract significant DEGs
#Here we are using a adjusted p-value of <= 0.05 and a log2FC of >= 1 or <= -1
#This log2FC is equivalent to a gene doubling or halving in expression
sig <- na.omit(resfull[resfull$padj <= 0.05 & resfull$log2FC >= 1 | 
  resfull$padj <= 0.05 & resfull$log2FC <= -1,])
#Count number of DEGs for each comparison
table(sig$sampleA,sig$sampleB)
#For each pairwise comparison, extract the DEGs
CvTsig <- sig[sig$sampleA == "Cabbage" & sig$sampleB == "TO1000",]
KvTsig <- sig[sig$sampleA == "Kale" & sig$sampleB == "TO1000",]
KvCsig <- sig[sig$sampleA == "Kale" & sig$sampleB == "Cabbage",]
#Output KvT DEGs
write.table(
  geneEx[geneEx$Gene %in% KvTsig$id,c(1,2,3,4,7,8,9,10,12,13,14)],
  paste(path,"KvT_DEGs.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE
)
#Output KvC DEGs
write.table(
  geneEx[geneEx$Gene %in% KvCsig$id,c(1,2,3,4,5,6,10,11,15,16)],
  paste(path,"KvC_DEGs.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE
)
#Output CvT DEGs
write.table(
  geneEx[geneEx$Gene %in% CvTsig$id,c(1,5,6,7,8,9,11,12,17,18)],
  paste(path,"CvT_DEGs.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE
)


#Create a table of DEG Numbers
DEGs <- data.frame(Comparison=c("Kale vs T01000","Kale vs Cabbage","Cabbage vs TO1000"),
  "Higher Expressed DEGs"=c(nrow(KvTsig[KvTsig$log2FC >= 1,]),
    nrow(KvCsig[KvCsig$log2FC >= 1,]),
    nrow(CvTsig[CvTsig$log2FC >= 1,])),
  "Lower Expressed DEGs"=c(nrow(KvTsig[KvTsig$log2FC <= -1,]),
    nrow(KvCsig[KvCsig$log2FC <= -1,]),
    nrow(CvTsig[CvTsig$log2FC <= -1,])),
  "Total DEGs"=c(nrow(KvTsig),nrow(KvCsig),nrow(CvTsig))
  )
colnames(DEGs) <- gsub("\\."," ",colnames(DEGs))
#Formatt and export the table with flextable
DEGsFT <- flextable(DEGs)
DEGsFT <- autofit(DEGsFT)
DEGsFT <- align(DEGsFT,align="center",part="all")
DEGsFT <- add_header_lines(DEGsFT,values=("Differentially Expressed Genes"))
DEGsFT <- bold(DEGsFT,part="header")
DEGsFT <- bold(DEGsFT,j="Comparison")
DEGsFT <- font(DEGsFT,fontname="Arial")
DEGsFT <- fontsize(DEGsFT,size=12,part="all")
save_as_docx(DEGsFT, path=paste(path,"DEGs_numbers.docx",sep=""))


#Make a venn diagram of DEGs
#First make a table of overlap between pairwise comparisons
d2 <- data.frame(id=unique(sig$id))
d2 <- data.frame(
	id=d2$id,
	KvT=ifelse(d2$id %in% sig[sig$sampleA == "Kale" & sig$sampleB == "TO1000",]$id,1,0),
	KvC=ifelse(d2$id %in% sig[sig$sampleA == "Kale" & sig$sampleB == "Cabbage",]$id,1,0),
	CvT=ifelse(d2$id %in% sig[sig$sampleA == "Cabbage" & sig$sampleB == "TO1000",]$id,1,0)
	)
#Make the venn diagram using VennDiagram and save as a pdf
library(VennDiagram)
cairo_pdf(paste(path,"DEG_Overlap.pdf",sep=""),width=6,height=6,family="Arial")
	draw.triple.venn(
		area1=nrow(subset(d2,KvT==1)),
		area2=nrow(subset(d2,KvC==1)),
		area3=nrow(subset(d2,CvT==1)),
		n12=nrow(subset(d2,KvT==1 & KvC==1)),
		n23=nrow(subset(d2,KvC==1 & CvT==1)),
		n13=nrow(subset(d2,KvT==1 & CvT==1)),
		n123=nrow(subset(d2,KvT==1 & KvC==1 & CvT==1)),
		category=c("Kale vs TO1000","Kale vs Cabbage","Cabbage vs TO1000"),
		lty="blank",fill=c("skyblue","pink1","mediumorchid"),
    fontfamily=rep("Arial"),cat.fontfamily=rep("Arial"),
    margin=0.08
	)
dev.off()


#Identify DEGs in common in both Kale comparisons
Kshared <- subset(d2,KvT==1 & KvC==1 & CvT==0)
KsharedUp <- geneEx[geneEx$Gene %in% Kshared$id & 
  geneEx$KvT_log2FC > 1 & geneEx$KvC_log2FC > 1,]
KsharedDown <- geneEx[geneEx$Gene %in% Kshared$id & 
  geneEx$KvT_log2FC < -1 & geneEx$KvC_log2FC < -1,]
#Read in descriptive annotations
BoAnnot <- read.csv("../misc/Bo_annotations.tsv",header=TRUE)
#Identify shared genes with increased or decreased expression, map these to
#the descriptive annotations, and output as csv files
KsharedUpAnnot <- merge(KsharedUp,BoAnnot)
write.csv(KsharedUpAnnot[c("Gene","Description","BLAST_hit")],
  paste(path,"Kale_shared_Up_Annotations.tsv",sep=""),quote=FALSE,row.names=FALSE)
KsharedDownAnnot <- merge(KsharedDown,BoAnnot)
write.csv(KsharedDownAnnot[c("Gene","Description","BLAST_hit")],
  paste(path,"Kale_shared_Down_Annotations.tsv",sep=""),quote=FALSE,row.names=FALSE)

#Compare Syntenic vs Non-Syntenic Genes between B. oleracea and Arabidopsis for 
#enrichment in DEGs
#Read in data and format tables
syn <- read.table("../misc/Bo-At-syntelogs.tsv",header=T,sep="\t")
synRes <- merge(resfull,syn,by.x="id",by.y="Bo_gene")
synSig <- merge(sig,syn,by.x="id",by.y="Bo_gene")

pSyn=data.frame(
  row.names=c("Syntenic: Genome","Syntenic: KvT DEGs","Syntenic: KvC DEGs",
    "Syntenic: CvT DEGs","Non-Syntenic: Genome","Non-Syntenic: KvT DEGs",
    "Non-Syntenic: KvC DEGs","Non-Syntenic: CvT DEGs"),
  percent=c(length(syn$Bo_gene %in% geneEx$Gene)/59225,
    nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="Kale" & sig$sampleB=="TO1000",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
    nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="Kale" & sig$sampleB=="Cabbage",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
    nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])/
      nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",]),
    (59225-length(syn$Bo_gene %in% geneEx$Gene))/59225,
    nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="Kale" & sig$sampleB=="TO1000",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
    nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="Kale" & sig$sampleB=="Cabbage",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
    nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])/
      nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])),
  order=c(1,2,3,4,5,6,7,8))
#Make Plot
ggsave(paste(path,"Syntenic_genes.pdf",sep=""),
  ggplot(pSyn) + 
    geom_bar(aes(x=reorder(row.names(pSyn),order),
     y=percent,fill=reorder(row.names(pSyn),order)),stat="identity") + 
    theme(panel.background=element_blank(),
      axis.line=element_line(color="black"),
      axis.text=element_text(size=12,color="black"),
      axis.title=element_text(size=18,color="black",face="bold"),
      axis.text.x=element_text(angle=270,hjust=0,vjust=0.5),
      plot.title = element_text(size=18,hjust = 0.5,color="black",face="bold"),
      legend.position="none") + 
    scale_y_continuous(expand=c(0,0),limits=c(0,0.7),labels=percent) + 
    scale_fill_manual(values=c("tomato2","dodgerblue3","palegreen4","khaki3",
      "tomato2","dodgerblue3","palegreen4","khaki3")) +
    ylab("Percentage of Genes") +
    xlab("") + ggtitle("Genes Syntenic to Arabidopsis")
)

#KvT subgenome
KvTsub=data.frame(
  row.names=c("LF Genome","LF DEGs","MF1 Genome","MF1 DEGs","MF2 Genome","MF2 DEGs"),
  percent=c(nrow(syn[syn$subgenome=="LF",])/59225,
    nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="LF",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
    nrow(syn[syn$subgenome=="MF1",])/59225,
    nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF1",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
    nrow(syn[syn$subgenome=="MF2",])/59225,
    nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF2",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",])),
  order=c(1,2,3,4,5,6))
#Make Plot
ggsave(paste(path,"KvT_subgenome.pdf",sep=""),
       ggplot(KvTsub) + 
         geom_bar(aes(x=reorder(row.names(KvTsub),order),
          y=percent,fill=reorder(row.names(KvTsub),order)),
          stat="identity") + 
         theme(panel.background=element_blank(),
          axis.line=element_line(color="black"),
          axis.text=element_text(size=12,color="black"),
          axis.title=element_text(size=18,color="black",face="bold"),
          axis.text.x=element_text(angle=270,hjust=0,vjust=0.5),
          plot.title = element_text(size=18,hjust = 0.5,color="black",face="bold"),
          legend.position="none") + 
         scale_y_continuous(expand=c(0,0),limits=c(0,0.23),labels=percent) + 
         scale_fill_manual(values=c("tomato2","tomato2","dodgerblue3",
          "dodgerblue3","black","black")) +
         ylab("Percentage of Genes") +
         xlab("") + ggtitle("Kale vs T01000 Subgenomes")
)
#KvC subgenome
KvCsub=data.frame(
  row.names=c("LF Genome","LF DEGs","MF1 Genome","MF1 DEGs","MF2 Genome","MF2 DEGs"),
  percent=c(nrow(syn[syn$subgenome=="LF",])/59225,
    nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="Cabbage" & synSig$subgenome=="LF",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
    nrow(syn[syn$subgenome=="MF1",])/59225,
    nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="Cabbage" & synSig$subgenome=="MF1",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
    nrow(syn[syn$subgenome=="MF2",])/59225,
    nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="Cabbage" & synSig$subgenome=="MF2",])/
      nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",])),
  order=c(1,2,3,4,5,6))
#Make Plot
ggsave(paste(path,"KvC_subgenome.pdf",sep=""),
       ggplot(KvCsub) + 
         geom_bar(aes(x=reorder(row.names(KvCsub),order),
          y=percent,fill=reorder(row.names(KvCsub),order)),
          stat="identity") + 
         theme(panel.background=element_blank(),
          axis.line=element_line(color="black"),
          axis.text=element_text(size=12,color="black"),
          axis.title=element_text(size=18,color="black",face="bold"),
          axis.text.x=element_text(angle=270,hjust=0,vjust=0.5),
          plot.title = element_text(size=18,hjust = 0.5,color="black",face="bold"),
          legend.position="none") + 
         scale_y_continuous(expand=c(0,0),limits=c(0,0.23),labels=percent) + 
         scale_fill_manual(values=c("tomato2","tomato2","dodgerblue3",
          "dodgerblue3","black","black")) +
         ylab("Percentage of Genes") +
         xlab("") + ggtitle("Kale vs Cabbage Subgenomes")
)
#CvT subgenome
CvTsub=data.frame(
  row.names=c("LF Genome","LF DEGs","MF1 Genome","MF1 DEGs","MF2 Genome","MF2 DEGs"),
  percent=c(nrow(syn[syn$subgenome=="LF",])/59225,
    nrow(synSig[synSig$sampleA=="Cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="LF",])/
      nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",]),
    nrow(syn[syn$subgenome=="MF1",])/59225,
    nrow(synSig[synSig$sampleA=="Cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF1",])/
      nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",]),
    nrow(syn[syn$subgenome=="MF2",])/59225,
    nrow(synSig[synSig$sampleA=="Cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF2",])/
      nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])),
  order=c(1,2,3,4,5,6))
#Make Plot
ggsave(paste(path,"CvT_subgenome.pdf",sep=""),
       ggplot(CvTsub) + 
         geom_bar(aes(x=reorder(row.names(CvTsub),order),
          y=percent,fill=reorder(row.names(CvTsub),order)),
          stat="identity") + 
         theme(panel.background=element_blank(),
          axis.line=element_line(color="black"),
          axis.text=element_text(size=12,color="black"),
          axis.title=element_text(size=18,color="black",face="bold"),
          axis.text.x=element_text(angle=270,hjust=0,vjust=0.5),
          plot.title = element_text(size=18,hjust = 0.5,color="black",face="bold"),
          legend.position="none") + 
         scale_y_continuous(expand=c(0,0),limits=c(0,0.23),labels=percent) + 
         scale_fill_manual(values=c("tomato2","tomato2","dodgerblue3",
          "dodgerblue3","black","black")) +
         ylab("Percentage of Genes") +
         xlab("") + ggtitle("Cabbage vs T01000 Subgenomes")
)


##GO term enrichment using topGO
#Create output directory
path <- "GO_terms/"
if(!dir.exists(path)){
  dir.create(path)
}
#Load libraries
library(topGO)
library(GO.db)
#Read in gene to GO term mapping in topGO format
goTerms <- readMappings(file="../misc/topGO.txt")
#Analyze GO terms for increased expression genes in Kale vs TO1000
#Format input table, selecting only genes with increasing expression
KvTgotermUP <- factor(as.integer(resKvT$id %in% KvTsig[KvTsig$log2FC > 1,]$id))
names(KvTgotermUP) <- resKvT$id
#Perform GO term enrichment
#topGO is a user defined function
#It combines multiple steps for each category of GO term and outputs a formatted table
#Here I am restricting analysis to only GO terms with ...
#at least 5 genes mapping to that to that term in the entire gene list
KvTgotermUP <- topGO(KvTgotermUP,goTerms,nodeSize=5,fdr=0.05,filename="KvT_up",
  path=path,writeData=TRUE)
#Analyze GO term enrichment for decreased expression genes in Kale vs TO1000
KvTgotermDOWN <- factor(as.integer(resKvT$id %in% KvTsig[KvTsig$log2FC < -1,]$id))
names(KvTgotermDOWN) <- resKvT$id
KvTgotermDOWN <- topGO(KvTgotermDOWN,goTerms,nodeSize=5,fdr=0.05,filename="KvT_down",
  path=path,writeData=TRUE)
#Analyze GO terms for increased expression genes in Kale vs Cabbage
KvCgotermUP <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC > 1,]$id))
names(KvCgotermUP) <- resKvC$id
KvCgotermUP <- topGO(KvCgotermUP,goTerms,nodeSize=5,fdr=0.05,filename="KvC_up",
  path=path,writeData=TRUE)
#Analyze GO term enrichment for decreased expression genes in Kale vs Cabbage
KvCgotermDOWN <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC < -1,]$id))
names(KvCgotermDOWN) <- resKvC$id
KvCgotermDOWN <- topGO(KvCgotermDOWN,goTerms,nodeSize=5,fdr=0.05,filename="KvC_down",
  path=path,writeData=TRUE)
#Analyze GO terms for increased expression genes in Cabbage vs TO1000 
CvTgotermUP <- factor(as.integer(resCvT$id %in% CvTsig[CvTsig$log2FC > 1,]$id))
names(CvTgotermUP) <- resCvT$id
CvTgotermUP <- topGO(CvTgotermUP,goTerms,nodeSize=5,fdr=0.05,filename="CvT_up",
  path=path,writeData=TRUE)
#Analyze GO term enrichment for decreased expression genes in Cabbage vs TO1000
CvTgotermDOWN <- factor(as.integer(resCvT$id %in% CvTsig[CvTsig$log2FC < -1,]$id))
names(CvTgotermDOWN) <- resCvT$id
CvTgotermDOWN <- topGO(CvTgotermDOWN,goTerms,nodeSize=5,fdr=0.05,filename="CvT_down",
  path=path,writeData=TRUE)
#Analyze GO terms for increased expression genes in Kale shared
KsharedGOtermUP <- factor(as.integer(geneEx$Gene %in% KsharedUp$Gene))
names(KsharedGOtermUP) <- geneEx$Gene
KsharedGOtermUP <- topGO(KsharedGOtermUP,goTerms,nodeSize=5,fdr=0.05,
  filename="Kale_shared_up",path=path,writeData=TRUE)
#Analyze GO term enrichment for decreased expression genes in Kale shared
KsharedGOtermDOWN <- factor(as.integer(geneEx$Gene %in% KsharedDown$Gene))
names(KsharedGOtermDOWN) <- geneEx$Gene
KsharedGOtermDOWN <- topGO(KsharedGOtermDOWN,goTerms,nodeSize=5,fdr=0.05,
  filename="Kale_shared_down",path=path,writeData=TRUE)


#KEGG analysis
#Check if directory "kegg" exists and if not create it
path <- "Kegg/"
if(!dir.exists(path)){
  dir.create(path)
}
#Read in mappings to Bo ncbi gene names
#These are the ones supported by KEGG for B. oleracea
ncbi <- read.csv("../misc/Bo2ncbi.csv",header=TRUE)
#Merge those mappings with Kale vs TO1000 results
KvTncbi <- merge(ncbi,resKvT,by.x="Boleracea_gene",by.y="id")
#Separate those results into genes with increased & decreased expresion
KvTncbiSigUp <- na.omit(KvTncbi[KvTncbi$padj < 0.05 & KvTncbi$log2FC >= 1,])
KvTncbiSigDown <- na.omit(KvTncbi[KvTncbi$padj < 0.05 & KvTncbi$log2FC <= -1,])
#Merge those mappings with Kale vs Cabbage results
KvCncbi <- merge(ncbi,resKvC,by.x="Boleracea_gene",by.y="id")
#Separate those results into genes with increased & decreased expresion
KvCncbiSigUp <- na.omit(KvCncbi[KvCncbi$padj < 0.05 & KvCncbi$log2FC >= 1,])
KvCncbiSigDown <- na.omit(KvCncbi[KvCncbi$padj < 0.05 & KvCncbi$log2FC <= -1,])
#Merge those mappings with Cabbage vs TO1000 results
CvTncbi <- merge(ncbi,resCvT,by.x="Boleracea_gene",by.y="id")
#Separate those results into genes with increased & decreased expresion
CvTncbiSigUp <- na.omit(CvTncbi[CvTncbi$padj < 0.05 & CvTncbi$log2FC >= 1,])
CvTncbiSigDown <- na.omit(CvTncbi[CvTncbi$padj < 0.05 & CvTncbi$log2FC <= -1,])
#Merge those mappings with Kale shared results for increased/decreased expression
KncbiSharedUp <- merge(ncbi,KsharedUp,by.x="Boleracea_gene",by.y="Gene")
KncbiSharedDown <- merge(ncbi,KsharedDown,by.x="Boleracea_gene",by.y="Gene")
#Use the package clusterProfiler to do KEGG enrichment
library(clusterProfiler)
#Perform enrichment test for Kale vs TO1000 genes with increased expression
#and export the results to table
KvTupKEGG <- enrichKEGG(KvTncbiSigUp$NCBI_gene, organism="boe")@result
write.table(KvTupKEGG[KvTupKEGG$p.adjust < 0.05,],
  paste(path,"KvT_up.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
#Perform enrichment test for Kale vs TO1000 genes with decreased expression
#and export the results to table
KvTdownKEGG <- enrichKEGG(KvTncbiSigDown$NCBI_gene, organism="boe")@result
write.table(KvTdownKEGG[KvTdownKEGG$p.adjust < 0.05,],
  paste(path,"KvT_down.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
#Perform enrichment test for Kale vs Cabbage genes with increased expression
#and export the results to table
KvCupKEGG <- enrichKEGG(KvCncbiSigUp$NCBI_gene, organism="boe")@result
write.table(KvCupKEGG[KvCupKEGG$p.adjust < 0.05,],
  paste(path,"KvC_up.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
#Perform enrichment test for Kale vs Cabbage genes with decreased expression
#and export the results to table
KvCdownKEGG <- enrichKEGG(KvCncbiSigDown$NCBI_gene, organism="boe")@result
write.table(KvCdownKEGG[KvCdownKEGG$p.adjust < 0.05,],
  paste(path,"KvC_down.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
#Perform enrichment test for Cabbage vs TO1000 genes with increased expression
#and export the results to table
CvTupKEGG <- enrichKEGG(CvTncbiSigUp$NCBI_gene, organism="boe")@result
write.table(CvTupKEGG[CvTupKEGG$p.adjust < 0.05,],
  paste(path,"CvT_up.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
#Perform enrichment test for Cabbage vs TO1000 genes with decreased expression
#and export the results to table
CvTdownKEGG <- enrichKEGG(CvTncbiSigDown$NCBI_gene, organism="boe")@result
write.table(CvTdownKEGG[CvTdownKEGG$p.adjust < 0.05,],
  paste(path,"CvT_down.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
#Perform enrichment test for Kale shared genes with increased expression
#and export the results to table
KsharedUpKEGG <- enrichKEGG(KncbiSharedUp$NCBI_gene, organism="boe")@result
write.table(KsharedUpKEGG[KsharedUpKEGG$p.adjust < 0.05,],
  paste(path,"Kale_shared_up.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
#Perform enrichment test for Kale shared genes with decreased expression
#and export the results to table
KsharedDownKEGG <- enrichKEGG(KncbiSharedDown$NCBI_gene, organism="boe")@result
write.table(KsharedDownKEGG[KsharedDownKEGG$p.adjust < 0.05,],
  paste(path,"Kale_shared_down.csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)


#Plot out data for specific genes of interest
#Create output directory
path <- "Genes_of_interest/"
if(!dir.exists(path)){
  dir.create(path)
}
#Read in gene list
goi <- read.csv("../misc/genes_of_interest.csv",header=TRUE)
#Iterate over each gene in that list
for(gene in goi$gene){
    #Use "tryCatch" to handle errors
    tryCatch({
        #Make a table of that gene's data from deseq2 using the plotCounts function
        x <- plotCounts(dds,gene,intgroup="condition",normalized=TRUE,transform=FALSE,
          returnData=TRUE)
        #Plot using ggplot2
        p <- ggplot(x) +
             geom_point(aes(x=condition,y= count,color=condition)) +
             theme(panel.background=element_blank(),
              axis.line=element_line(color="black"),
              axis.text=element_text(color="black"),
              axis.title=element_text(color="black",face="bold"),
              legend.position="none") + xlab("Genotype")
    #Save as a pdf
    ggsave(paste(path,gene,"_counts.pdf",sep=""),p,width=5,height=4)
    #How to handle potential errors
    },error=function(e){})
}
