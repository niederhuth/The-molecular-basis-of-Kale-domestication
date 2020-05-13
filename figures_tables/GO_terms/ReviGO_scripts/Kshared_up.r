

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
revigo.data <- rbind(c("GO:0009657","plastid organization",0.024,1.4857,0.832,0.000,"plastid organization"),
c("GO:0043271","negative regulation of ion transport",0.022,1.4236,0.729,0.000,"negative regulation of ion transport"),
c("GO:0050898","nitrile metabolic process",0.000,2.0543,0.803,0.000,"nitrile metabolism"),
c("GO:0080028","nitrile biosynthetic process",0.000,2.0891,0.774,0.103,"nitrile metabolism"),
c("GO:0051707","response to other organism",0.299,1.4857,0.765,0.000,"response to other organism"),
c("GO:0010106","cellular response to iron ion starvation",0.002,1.4857,0.758,0.506,"response to other organism"),
c("GO:0008299","isoprenoid biosynthetic process",0.442,1.6174,0.633,0.024,"isoprenoid biosynthesis"),
c("GO:0008610","lipid biosynthetic process",2.123,1.4169,0.627,0.696,"isoprenoid biosynthesis"),
c("GO:1901570","fatty acid derivative biosynthetic process",0.009,1.4857,0.709,0.268,"isoprenoid biosynthesis"),
c("GO:0006629","lipid metabolic process",3.522,1.4236,0.732,0.179,"isoprenoid biosynthesis"),
c("GO:1901568","fatty acid derivative metabolic process",0.017,1.4857,0.762,0.111,"isoprenoid biosynthesis"),
c("GO:0000272","polysaccharide catabolic process",0.288,1.4236,0.799,0.051,"polysaccharide catabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="GO_terms/ReviGO_plots/Kshared_up_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the treemap command documentation for all possible parameters - there are a lot more
treemap(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "Kale Shared Higher Expressed Genes",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
