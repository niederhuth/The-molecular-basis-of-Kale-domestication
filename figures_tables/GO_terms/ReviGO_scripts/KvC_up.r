

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
revigo.data <- rbind(c("GO:0018874","benzoate metabolic process",0.004,2.2055,0.804,0.000,"benzoate metabolism"),
c("GO:0043271","negative regulation of ion transport",0.022,1.4758,0.829,0.000,"negative regulation of ion transport"),
c("GO:0080027","response to herbivore",0.000,2.9987,0.671,0.000,"response to herbivore"),
c("GO:0006950","response to stress",4.575,1.4873,0.660,0.501,"response to herbivore"),
c("GO:0034052","positive regulation of plant-type hypersensitive response",0.000,2.2055,0.608,0.141,"response to herbivore"),
c("GO:0009605","response to external stimulus",1.370,1.6879,0.666,0.220,"response to herbivore"),
c("GO:0016046","detection of fungus",0.000,2.2055,0.685,0.422,"response to herbivore"),
c("GO:0010106","cellular response to iron ion starvation",0.002,1.6665,0.636,0.362,"response to herbivore"),
c("GO:0007032","endosome organization",0.020,1.5815,0.849,0.017,"endosome organization"),
c("GO:0050898","nitrile metabolic process",0.000,2.1650,0.840,0.019,"nitrile metabolism"),
c("GO:0080028","nitrile biosynthetic process",0.000,2.2055,0.830,0.103,"nitrile metabolism"),
c("GO:0015979","photosynthesis",0.183,2.2055,0.841,0.036,"photosynthesis"),
c("GO:1901568","fatty acid derivative metabolic process",0.017,2.1650,0.822,0.083,"fatty acid derivative metabolism"),
c("GO:0015977","carbon fixation",0.036,1.4873,0.821,0.095,"carbon fixation"));

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
