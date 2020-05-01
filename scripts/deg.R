#Import user-defined functions
source("../scripts/functions.R")

if(!dir.exists("paper")){
	dir.create("paper")
}

#Make pdf table of mapping statistics
library(scales)
library(flextable)
mapStats <- read.csv("Mapping_Stats.csv",header=TRUE)
#Add a percentage mapped column
mapStats$Pecent.Uniquely.Mapped <- percent(mapStats$Uniquely.Mapped/mapStats$Trimmed.Reads,
											accuracy=0.1)
#Get rid of period in the column names
colnames(mapStats) <- gsub("\\."," ",colnames(mapStats))
#Create table with flextable
myft1 <- flextable(mapStats)
myft1 <- autofit(myft1)
myft1 <- align(myft1,align="center",part="all")
myft1 <- add_header_lines(myft1,values=("Supplementary Table 1: Mapping Statistics"))
myft1 <- bold(myft1,part="header")
myft1 <- bold(myft1,j="Sample")
myft1 <- font(myft1,fontname="Arial")
myft1 <- fontsize(myft1,size=12,part="all")
save_as_docx(myft1, path="paper/Supplementary_Table_1.docx")


#Run DESeq2 analysis
#Load Deseq2
library(DESeq2)
#Read in sample table
sampleTable <- read.csv("../misc/deseq2_samples.csv",header=T)
#Read in counts table
dds <- DESeqDataSetFromHTSeqCount(sampleTable,design= ~ condition)
#Set reference level
dds$condition <- relevel(dds$condition, ref="TO1000")
#Estimate size factors
dds <- estimateSizeFactors(dds)
#Estimate dispersions
dds <- estimateDispersions(dds,fitType="parametric")
#Run Wald test
dds <- nbinomWaldTest(dds)

#Make diagnostic figures of samples
#Make a heat map of samples
rld <- rlog(dds, blind=TRUE)
pdf("Sample_Comparison_Heat_Map.pdf",width=6,height=6,paper='special')
	#sampleHeatMap is a user-defined function
	sampleHeatMap(rld)
dev.off()
#Make a PCA plot of samples
pdf("Sample_Comparison_PCA.pdf",width=6,height=6,paper='special')
	#pcaPlot is a user-defined function
	pcaPlot(rld)
dev.off()

#Make results tables for each pairwise comparison
#makeResultsTable is a user-defined function
resKvT <- makeResultsTable(dds,"Kale","TO1000",lfcThreshold=0,filter=F)
resKvC <- makeResultsTable(dds,"Kale","Cabbage",lfcThreshold=0,filter=F)
resCvT <- makeResultsTable(dds,"Cabbage","TO1000",lfcThreshold=0,filter=F)
#Combine results tables
resfull <- as.data.frame(rbind(resKvT,resKvC,resCvT))
#Adjust p-values for all results
#This is necessary because the adjusted p-value is only for the individual pairwise comparisons
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
write.table(geneEx,"Gene_Expression_Table.tsv",sep="\t",quote=FALSE,row.names=FALSE)

#Extract significant DEGs
#Here we are using a adjusted p-value of <= 0.05 and a log2FC of >= 1 or <= -1
#This log2FC is equivalent to a gene doubling or halving in expression
sig <- na.omit(resfull[resfull$padj <= 0.05 & resfull$log2FC >= 1 | resfull$padj <= 0.05 & resfull$log2FC <= -1,])
#Count number of DEGs for each comparison
table(sig$sampleA,sig$sampleB)
#For each pairwise comparison, extract the DEGs
CvTsig <- sig[sig$sampleA == "Cabbage" & sig$sampleB == "TO1000",]
KvTsig <- sig[sig$sampleA == "Kale" & sig$sampleB == "TO1000",]
KvCsig <- sig[sig$sampleA == "Kale" & sig$sampleB == "Cabbage",]

#Create a table of DEG Numbers
df3 <- data.frame(Comparison=c("Kale vs T01000","Kale vs Cabbage","Cabbage vs TO1000"),
  "Higher Expressed DEGs"=c(nrow(KvTsig[KvTsig$log2FC >= 1,]),
    nrow(KvCsig[KvCsig$log2FC >= 1,]),
    nrow(CvTsig[CvTsig$log2FC >= 1,])),
  "Lower Expressed DEGs"=c(nrow(KvTsig[KvTsig$log2FC <= -1,]),
    nrow(KvCsig[KvCsig$log2FC <= -1,]),
    nrow(CvTsig[CvTsig$log2FC <= -1,])),
  "Total DEGs"=c(nrow(KvTsig),nrow(KvCsig),nrow(CvTsig))
  )
colnames(df3) <- gsub("\\."," ",colnames(df3))
#Formatt and export the table with flextable
myft2 <- flextable(df3)
myft2 <- autofit(myft2)
myft2 <- align(myft2,align="center",part="all")
myft2 <- add_header_lines(myft2,values=("A: Differentially Expressed Genes"))
myft2 <- bold(myft2,part="header")
myft2 <- bold(myft2,j="Comparison")
myft2 <- font(myft2,fontname="Arial")
myft2 <- fontsize(myft2,size=12,part="all")
save_as_docx(myft2, path="paper/Table_1.docx")

#Make a venn diagram of DEGs
#First make a table of overlap between pairwise comparisons
d2 <- data.frame(id=unique(sig$id))
d2 <- data.frame(
	id=d2$id,
	KvT=ifelse(d2$id %in% sig[sig$sampleA == "Kale" & sig$sampleB == "TO1000",]$id, 1, 0),
	KvC=ifelse(d2$id %in% sig[sig$sampleA == "Kale" & sig$sampleB == "Cabbage",]$id, 1, 0),
	CvT=ifelse(d2$id %in% sig[sig$sampleA == "Cabbage" & sig$sampleB == "TO1000",]$id, 1, 0)
	)
#Make the venn diagram using VennDiagram and save as a pdf
library(VennDiagram)
pdf("DEG_Overlap.pdf",width=6,height=6,paper='special')
	draw.triple.venn(
		area1=nrow(subset(d2,CvT==1)),
		area2=nrow(subset(d2,KvT==1)),
		area3=nrow(subset(d2,KvC==1)),
		n12=nrow(subset(d2,CvT==1 & KvT==1)),
		n23=nrow(subset(d2,KvT==1 & KvC==1)),
		n13=nrow(subset(d2,CvT==1 & KvC==1)),
		n123=nrow(subset(d2,CvT==1 & KvT==1 & KvC==1)),
		category=c("Cabbage v TO1000","Kale v TO100","Kale v Cabbage"),
			lty="blank",fill=c("skyblue","pink1","mediumorchid")
	)
dev.off()

##GO term enrichment using topGO
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
KvTgotermUP <- topGO(KvTgotermUP,goTerms,nodeSize=5,"KvT_up",writeData=TRUE)
#Analyze GO term enrichment for decreased expression genes in Kale vs TO1000
KvTgotermDOWN <- factor(as.integer(resKvT$id %in% KvTsig[KvTsig$log2FC < -1,]$id))
names(KvTgotermDOWN) <- resKvT$id
KvTgotermDOWN <- topGO(KvTgotermDOWN,goTerms,nodeSize=5,"KvT_down",writeData=TRUE)

#Analyze GO terms for increased expression genes in Kale vs Cabbage
KvCgotermUP <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC > 1,]$id))
names(KvCgotermUP) <- resKvC$id
KvCgotermUP <- topGO(KvCgotermUP,goTerms,nodeSize=5,"KvC_up",writeData=TRUE)
#Analyze GO term enrichment for decreased expression genes in Kale vs Cabbage
KvCgotermDOWN <- factor(as.integer(resKvC$id %in% KvCsig[KvCsig$log2FC < -1,]$id))
names(KvCgotermDOWN) <- resKvC$id
KvCgotermDOWN <- topGO(KvCgotermDOWN,goTerms,nodeSize=5,"KvC_down",writeData=TRUE)

#Analyze GO terms for increased expression genes in Cabbage vs TO1000 
CvTgotermUP <- factor(as.integer(resCvT$id %in% CvTsig[CvTsig$log2FC > 1,]$id))
names(CvTgotermUP) <- resCvT$id
CvTgotermUP <- topGO(CvTgotermUP,goTerms,nodeSize=5,"CvT_up",writeData=TRUE)
#Analyze GO term enrichment for decreased expression genes in Cabbage vs TO1000
CvTgotermDOWN <- factor(as.integer(resCvT$id %in% CvTsig[CvTsig$log2FC < -1,]$id))
names(CvTgotermDOWN) <- resCvT$id
CvTgotermDOWN <- topGO(CvTgotermDOWN,goTerms,nodeSize=5,"CvT_down",writeData=TRUE)


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
#library(clusterProfiler)
#xx <- enrichMKEGG(KvTsig$id, organism='boe', minGSSize=1)


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

#Compare Syntenic vs Non-Syntenic Genes between B. oleracea and Arabidopsis for 
#enrichment in DEGs
#Read in data and format tables
syn <- read.table("../misc/Bo-At-syntelogs.tsv",header=T,sep="\t")
synRes <- merge(resfull,syn,by.x="id",by.y="Bo_gene")
synSig <- merge(sig,syn,by.x="id",by.y="Bo_gene")

library(scales)
pSyn=data.frame(
  row.names=c("Genome: Percent Syntenic",
              "KvT DEGs: Percent Syntenic",
              "KvC DEGs: Percent Syntenic",
              "CvT DEGs: Percent Syntenic",
              "Genome: Percent Non-Syntenic",
              "KvT DEGs: Percent Non-Syntenic",
              "KvC DEGs: Percent Non-Syntenic",
              "CvT DEGs: Percent Non-Syntenic"),
  percent=c(length(syn$Bo_gene %in% geneEx$Gene)/nrow(geneEx),
            nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="Kale" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
            nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="Kale" & sig$sampleB=="Cabbage",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
            nrow(sig[sig$id %in% syn$Bo_gene & sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",]),
            (nrow(geneEx)-length(syn$Bo_gene %in% geneEx$Gene))/nrow(geneEx),
            nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="Kale" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
            nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="Kale" & sig$sampleB=="Cabbage",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
            nrow(sig[!(sig$id %in% syn$Bo_gene) & sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])/nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])),
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
  percent=c(nrow(syn[syn$subgenome=="LF",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="LF",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF1",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF1",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF2",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF2",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="TO1000",])),
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
tmp <- synRes[synRes$sampleA=="Kale" & synRes$sampleB=="TO1000",]
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
  percent=c(nrow(syn[syn$subgenome=="LF",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="Cabbage" & synSig$subgenome=="LF",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
            nrow(syn[syn$subgenome=="MF1",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="Cabbage" & synSig$subgenome=="MF1",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",]),
            nrow(syn[syn$subgenome=="MF2",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Kale" & synSig$sampleB=="Cabbage" & synSig$subgenome=="MF2",])/nrow(sig[sig$sampleA=="Kale" & sig$sampleB=="Cabbage",])),
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
tmp <- synRes[synRes$sampleA=="Kale" & synRes$sampleB=="Cabbage",]
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
  percent=c(nrow(syn[syn$subgenome=="LF",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="LF",])/nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF1",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF1",])/nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",]),
            nrow(syn[syn$subgenome=="MF2",])/nrow(geneEx),
            nrow(synSig[synSig$sampleA=="Cabbage" & synSig$sampleB=="TO1000" & synSig$subgenome=="MF2",])/nrow(sig[sig$sampleA=="Cabbage" & sig$sampleB=="TO1000",])),
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
tmp <- synRes[synRes$sampleA=="Cabbage" & synRes$sampleB=="TO1000",]
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
#Make Venn Diagram
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
