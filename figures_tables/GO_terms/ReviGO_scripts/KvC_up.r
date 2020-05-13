

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0050898","nitrile metabolic process",0.000,2.7076,0.875,0.000,"nitrile metabolism"),
c("GO:0080028","nitrile biosynthetic process",0.000,2.7704,0.863,0.103,"nitrile metabolism"),
c("GO:0080027","response to herbivore",0.000,2.7704,0.819,0.000,"response to herbivore"),
c("GO:0006950","response to stress",4.575,2.4693,0.818,0.238,"response to herbivore"),
c("GO:0051707","response to other organism",0.299,1.3724,0.793,0.634,"response to herbivore"),
c("GO:0009605","response to external stimulus",1.370,2.2652,0.820,0.501,"response to herbivore"),
c("GO:0015979","photosynthesis",0.183,2.4636,0.889,0.007,"photosynthesis"),
c("GO:0015977","carbon fixation",0.036,2.4413,0.866,0.021,"carbon fixation"),
c("GO:0008299","isoprenoid biosynthetic process",0.442,1.4113,0.839,0.118,"carbon fixation"),
c("GO:0007033","vacuole organization",0.102,1.3444,0.902,0.025,"vacuole organization"),
c("GO:0010629","negative regulation of gene expression",0.784,2.0434,0.623,0.035,"negative regulation of gene expression"),
c("GO:0032879","regulation of localization",0.726,1.5723,0.787,0.271,"negative regulation of gene expression"),
c("GO:0043271","negative regulation of ion transport",0.022,1.3444,0.625,0.579,"negative regulation of gene expression"),
c("GO:0010361","regulation of anion channel activity by blue light",0.000,1.3444,0.697,0.399,"negative regulation of gene expression"),
c("GO:0071902","positive regulation of protein serine/threonine kinase activity",0.075,1.5796,0.734,0.314,"negative regulation of gene expression"),
c("GO:0005975","carbohydrate metabolic process",5.260,1.6040,0.874,0.058,"carbohydrate metabolism"),
c("GO:0010021","amylopectin biosynthetic process",0.001,1.3724,0.835,0.068,"amylopectin biosynthesis"),
c("GO:0019759","glycosinolate catabolic process",0.001,1.3444,0.819,0.219,"amylopectin biosynthesis"),
c("GO:2000896","amylopectin metabolic process",0.001,1.3444,0.845,0.221,"amylopectin biosynthesis"),
c("GO:1901568","fatty acid derivative metabolic process",0.017,1.4113,0.867,0.095,"fatty acid derivative metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="GO_terms/ReviGO_plots/KvC_up_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the treemap command documentation for all possible parameters - there are a lot more
treemap(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "Kale vs Cabbage Higher Expressed Genes",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
