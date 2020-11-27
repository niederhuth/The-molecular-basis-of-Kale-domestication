

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
revigo.data <- rbind(c("GO:0044419","interspecies interaction between organisms",0.262,3.2108,0.907,0.000,"interspecies interaction between organisms"),
c("GO:0050896","response to stimulus",12.210,1.7213,0.918,0.000,"response to stimulus"),
c("GO:0050898","nitrile metabolic process",0.000,3.1798,0.877,0.000,"nitrile metabolism"),
c("GO:0080028","nitrile biosynthetic process",0.000,3.2108,0.853,0.103,"nitrile metabolism"),
c("GO:0080027","response to herbivore",0.000,2.8405,0.795,0.000,"response to herbivore"),
c("GO:0034052","positive regulation of plant-type hypersensitive response",0.000,1.6677,0.718,0.141,"response to herbivore"),
c("GO:0009607","response to biotic stimulus",0.342,1.5172,0.818,0.204,"response to herbivore"),
c("GO:0016046","detection of fungus",0.000,1.4592,0.802,0.422,"response to herbivore"),
c("GO:0010106","cellular response to iron ion starvation",0.002,1.5059,0.766,0.362,"response to herbivore"),
c("GO:0015979","photosynthesis",0.183,2.4661,0.871,0.007,"photosynthesis"),
c("GO:1901568","fatty acid derivative metabolic process",0.017,2.6422,0.813,0.020,"fatty acid derivative metabolism"),
c("GO:0008610","lipid biosynthetic process",2.123,2.1193,0.705,0.218,"fatty acid derivative metabolism"),
c("GO:0008299","isoprenoid biosynthetic process",0.442,1.5172,0.699,0.696,"fatty acid derivative metabolism"),
c("GO:0006629","lipid metabolic process",3.522,1.9201,0.762,0.131,"fatty acid derivative metabolism"),
c("GO:0006082","organic acid metabolic process",9.086,2.9987,0.712,0.366,"fatty acid derivative metabolism"),
c("GO:0044281","small molecule metabolic process",15.138,3.1798,0.761,0.300,"fatty acid derivative metabolism"),
c("GO:1901570","fatty acid derivative biosynthetic process",0.009,1.7147,0.783,0.300,"fatty acid derivative metabolism"),
c("GO:1901136","carbohydrate derivative catabolic process",0.423,1.5172,0.870,0.032,"carbohydrate derivative catabolism"),
c("GO:0006468","protein phosphorylation",4.137,1.5216,0.843,0.063,"protein phosphorylation"),
c("GO:0018874","benzoate metabolic process",0.004,1.9987,0.791,0.083,"benzoate metabolism"),
c("GO:0010119","regulation of stomatal movement",0.004,1.3033,0.824,0.092,"regulation of stomatal movement"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()


stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="GO_terms/ReviGO_plots/Kshared_up_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
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
