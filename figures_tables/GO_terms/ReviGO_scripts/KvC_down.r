

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
revigo.data <- rbind(c("GO:0010200","response to chitin",0.004,5.0000,0.721,0.000,"response to chitin"),
c("GO:0080167","response to karrikin",0.006,2.6210,0.671,0.187,"response to chitin"),
c("GO:0010117","photoprotection",0.000,1.4380,0.678,0.461,"response to chitin"),
c("GO:0010438","cellular response to sulfur starvation",0.000,3.0981,0.716,0.158,"response to chitin"),
c("GO:0019757","glycosinolate metabolic process",0.003,2.9670,0.606,0.000,"glycosinolate metabolism"),
c("GO:0045839","negative regulation of mitotic nuclear division",0.051,1.8852,0.685,0.104,"glycosinolate metabolism"),
c("GO:0060255","regulation of macromolecule metabolic process",11.716,1.7272,0.634,0.493,"glycosinolate metabolism"),
c("GO:0001932","regulation of protein phosphorylation",0.430,1.4380,0.639,0.208,"glycosinolate metabolism"),
c("GO:0010467","gene expression",19.671,2.3499,0.747,0.038,"gene expression"),
c("GO:0006790","sulfur compound metabolic process",1.822,2.2619,0.771,0.042,"sulfur compound metabolism"));

stuff <- data.frame(revigo.data);
stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="GO_terms/ReviGO_plots/KvC_down_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the treemap command documentation for all possible parameters - there are a lot more
treemap(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "Kale vs Cabbage Lower Expressed Genes",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
