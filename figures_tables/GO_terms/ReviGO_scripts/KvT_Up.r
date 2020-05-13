

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
revigo.data <- rbind(c("GO:0010035","response to inorganic substance",0.317,4.2062,0.828,0.000,"response to inorganic substance"),
c("GO:0009625","response to insect",0.001,1.7061,0.873,0.658,"response to inorganic substance"),
c("GO:0080184","response to phenylpropanoid",0.000,1.6540,0.868,0.369,"response to inorganic substance"),
c("GO:0001101","response to acid chemical",0.124,1.6832,0.836,0.542,"response to inorganic substance"),
c("GO:0009734","auxin-activated signaling pathway",0.035,1.4863,0.726,0.492,"response to inorganic substance"),
c("GO:0051707","response to other organism",0.299,3.7762,0.841,0.319,"response to inorganic substance"),
c("GO:0009611","response to wounding",0.127,1.9305,0.854,0.297,"response to inorganic substance"),
c("GO:0006970","response to osmotic stress",0.082,1.7372,0.837,0.669,"response to inorganic substance"),
c("GO:0009415","response to water",0.026,2.2443,0.817,0.482,"response to inorganic substance"),
c("GO:0006950","response to stress",4.575,2.5107,0.845,0.418,"response to inorganic substance"),
c("GO:0010363","regulation of plant-type hypersensitive response",0.001,1.7893,0.770,0.318,"response to inorganic substance"),
c("GO:0080027","response to herbivore",0.000,3.4671,0.876,0.634,"response to inorganic substance"),
c("GO:0015979","photosynthesis",0.183,8.3028,0.933,0.000,"photosynthesis"),
c("GO:0050896","response to stimulus",12.210,6.8665,0.973,0.000,"response to stimulus"),
c("GO:1901568","fatty acid derivative metabolic process",0.017,6.8665,0.884,0.009,"fatty acid derivative metabolism"),
c("GO:0044281","small molecule metabolic process",15.138,6.8665,0.829,0.300,"fatty acid derivative metabolism"),
c("GO:0019748","secondary metabolic process",0.138,2.5020,0.877,0.158,"fatty acid derivative metabolism"),
c("GO:0006629","lipid metabolic process",3.522,6.5114,0.837,0.131,"fatty acid derivative metabolism"),
c("GO:0006082","organic acid metabolic process",9.086,6.0809,0.779,0.366,"fatty acid derivative metabolism"),
c("GO:0015977","carbon fixation",0.036,2.7423,0.880,0.139,"fatty acid derivative metabolism"),
c("GO:0009657","plastid organization",0.024,1.7372,0.943,0.022,"plastid organization"),
c("GO:0042545","cell wall modification",0.045,1.4836,0.943,0.304,"plastid organization"),
c("GO:0080028","nitrile biosynthetic process",0.000,3.5864,0.913,0.029,"nitrile biosynthesis"),
c("GO:0016143","S-glycoside metabolic process",0.003,1.8884,0.827,0.381,"nitrile biosynthesis"),
c("GO:0050898","nitrile metabolic process",0.000,3.4413,0.933,0.120,"nitrile biosynthesis"),
c("GO:1901136","carbohydrate derivative catabolic process",0.423,2.0581,0.876,0.328,"nitrile biosynthesis"),
c("GO:0019682","glyceraldehyde-3-phosphate metabolic process",0.444,1.5097,0.788,0.579,"nitrile biosynthesis"),
c("GO:0019757","glycosinolate metabolic process",0.003,3.4413,0.746,0.115,"nitrile biosynthesis"),
c("GO:0055129","L-proline biosynthetic process",0.085,1.3184,0.791,0.252,"nitrile biosynthesis"),
c("GO:1901657","glycosyl compound metabolic process",2.961,1.6132,0.788,0.576,"nitrile biosynthesis"),
c("GO:1901615","organic hydroxy compound metabolic process",0.831,4.0680,0.937,0.033,"organic hydroxy compound metabolism"),
c("GO:0009765","photosynthesis, light harvesting",0.019,2.6843,0.918,0.040,"photosynthesis, light harvesting"),
c("GO:0010119","regulation of stomatal movement",0.004,3.0861,0.837,0.043,"regulation of stomatal movement"),
c("GO:0051252","regulation of RNA metabolic process",10.029,2.2443,0.765,0.443,"regulation of stomatal movement"),
c("GO:1903792","negative regulation of anion transport",0.004,2.6843,0.767,0.174,"regulation of stomatal movement"),
c("GO:0006639","acylglycerol metabolic process",0.041,1.7423,0.822,0.566,"regulation of stomatal movement"),
c("GO:0019725","cellular homeostasis",1.253,2.7076,0.745,0.173,"regulation of stomatal movement"),
c("GO:0044255","cellular lipid metabolic process",2.704,2.4693,0.767,0.234,"regulation of stomatal movement"),
c("GO:0065008","regulation of biological quality",3.395,2.1440,0.855,0.325,"regulation of stomatal movement"),
c("GO:0018298","protein-chromophore linkage",0.095,5.0132,0.922,0.045,"protein-chromophore linkage"),
c("GO:0007568","aging",0.088,1.7893,0.895,0.052,"aging"),
c("GO:0090693","plant organ senescence",0.006,1.8884,0.893,0.418,"aging"),
c("GO:0006790","sulfur compound metabolic process",1.822,1.7061,0.926,0.058,"sulfur compound metabolism"),
c("GO:1901570","fatty acid derivative biosynthetic process",0.009,3.5151,0.856,0.087,"fatty acid derivative biosynthesis"),
c("GO:1901617","organic hydroxy compound biosynthetic process",0.383,1.5097,0.912,0.111,"fatty acid derivative biosynthesis"),
c("GO:0010118","stomatal movement",0.006,2.0576,0.890,0.094,"stomatal movement"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="GO_terms/ReviGO_plots/KvT_up_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the treemap command documentation for all possible parameters - there are a lot more
treemap(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "Kale vs TO1000 Higher Expressed Genes",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
