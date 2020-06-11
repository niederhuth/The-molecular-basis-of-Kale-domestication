

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
revigo.data <- rbind(c("GO:0001906","cell killing",0.032,1.3444,0.932,0.000,"cell killing"),
c("GO:0043271","negative regulation of ion transport",0.022,1.6284,0.882,0.000,"negative regulation of ion transport"),
c("GO:0050896","response to stimulus",12.210,1.6868,0.940,0.000,"response to stimulus"),
c("GO:0050898","nitrile metabolic process",0.000,5.3028,0.891,0.000,"nitrile metabolism"),
c("GO:0080028","nitrile biosynthetic process",0.000,5.3665,0.876,0.103,"nitrile metabolism"),
c("GO:0080027","response to herbivore",0.000,3.5847,0.830,0.000,"response to herbivore"),
c("GO:0080184","response to phenylpropanoid",0.000,1.5363,0.870,0.148,"response to herbivore"),
c("GO:0006950","response to stress",4.575,1.6868,0.833,0.238,"response to herbivore"),
c("GO:0010106","cellular response to iron ion starvation",0.002,2.3096,0.825,0.362,"response to herbivore"),
c("GO:0051707","response to other organism",0.299,2.4652,0.799,0.634,"response to herbivore"),
c("GO:0009607","response to biotic stimulus",0.342,1.6868,0.843,0.421,"response to herbivore"),
c("GO:0015979","photosynthesis",0.183,2.0134,0.902,0.007,"photosynthesis"),
c("GO:1901568","fatty acid derivative metabolic process",0.017,3.5910,0.841,0.020,"fatty acid derivative metabolism"),
c("GO:0008610","lipid biosynthetic process",2.123,1.6132,0.753,0.696,"fatty acid derivative metabolism"),
c("GO:0019748","secondary metabolic process",0.138,2.0966,0.832,0.131,"fatty acid derivative metabolism"),
c("GO:1901570","fatty acid derivative biosynthetic process",0.009,2.4468,0.820,0.268,"fatty acid derivative metabolism"),
c("GO:0008299","isoprenoid biosynthetic process",0.442,2.4468,0.757,0.111,"fatty acid derivative metabolism"),
c("GO:0006629","lipid metabolic process",3.522,1.6868,0.794,0.179,"fatty acid derivative metabolism"),
c("GO:0006082","organic acid metabolic process",9.086,1.8215,0.755,0.366,"fatty acid derivative metabolism"),
c("GO:0044281","small molecule metabolic process",15.138,1.7687,0.787,0.300,"fatty acid derivative metabolism"),
c("GO:0015977","carbon fixation",0.036,1.4266,0.836,0.118,"fatty acid derivative metabolism"),
c("GO:0009657","plastid organization",0.024,1.7671,0.920,0.022,"plastid organization"),
c("GO:1901136","carbohydrate derivative catabolic process",0.423,2.6754,0.753,0.032,"carbohydrate derivative catabolism"),
c("GO:0019759","glycosinolate catabolic process",0.001,2.8673,0.674,0.338,"carbohydrate derivative catabolism"),
c("GO:0016143","S-glycoside metabolic process",0.003,1.6454,0.742,0.353,"carbohydrate derivative catabolism"),
c("GO:0044273","sulfur compound catabolic process",0.117,2.4126,0.743,0.473,"carbohydrate derivative catabolism"),
c("GO:1901658","glycosyl compound catabolic process",0.085,1.5453,0.689,0.573,"carbohydrate derivative catabolism"),
c("GO:0000272","polysaccharide catabolic process",0.288,1.6647,0.790,0.519,"carbohydrate derivative catabolism"));

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
