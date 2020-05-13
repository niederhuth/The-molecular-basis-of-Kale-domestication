

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
revigo.data <- rbind(c("GO:0006790","sulfur compound metabolic process",1.822,6.3595,0.904,0.000,"sulfur compound metabolism"),
c("GO:0009611","response to wounding",0.127,4.2652,0.801,0.000,"response to wounding"),
c("GO:0009625","response to insect",0.001,2.8278,0.832,0.337,"response to wounding"),
c("GO:1901700","response to oxygen-containing compound",0.503,1.6868,0.773,0.611,"response to wounding"),
c("GO:0010438","cellular response to sulfur starvation",0.000,3.0714,0.809,0.294,"response to wounding"),
c("GO:0034599","cellular response to oxidative stress",0.224,1.3221,0.738,0.569,"response to wounding"),
c("GO:0051707","response to other organism",0.299,2.0714,0.781,0.658,"response to wounding"),
c("GO:0006979","response to oxidative stress",0.575,2.9999,0.783,0.498,"response to wounding"),
c("GO:0006970","response to osmotic stress",0.082,4.0434,0.773,0.427,"response to wounding"),
c("GO:0010200","response to chitin",0.004,1.9050,0.809,0.428,"response to wounding"),
c("GO:0009414","response to water deprivation",0.022,3.2601,0.753,0.662,"response to wounding"),
c("GO:0010117","photoprotection",0.000,1.3221,0.834,0.523,"response to wounding"),
c("GO:0006950","response to stress",4.575,2.4508,0.783,0.379,"response to wounding"),
c("GO:0042221","response to chemical",3.071,1.4352,0.788,0.562,"response to wounding"),
c("GO:0010035","response to inorganic substance",0.317,2.1545,0.779,0.418,"response to wounding"),
c("GO:0050896","response to stimulus",12.210,6.5114,0.969,0.000,"response to stimulus"),
c("GO:0044550","secondary metabolite biosynthetic process",0.101,6.0168,0.720,0.013,"secondary metabolite biosynthesis"),
c("GO:1901576","organic substance biosynthetic process",30.365,2.9999,0.807,0.282,"secondary metabolite biosynthesis"),
c("GO:0032774","RNA biosynthetic process",10.925,3.0185,0.739,0.462,"secondary metabolite biosynthesis"),
c("GO:1901617","organic hydroxy compound biosynthetic process",0.383,1.9050,0.865,0.158,"secondary metabolite biosynthesis"),
c("GO:0044283","small molecule biosynthetic process",5.677,1.5079,0.710,0.655,"secondary metabolite biosynthesis"),
c("GO:0051252","regulation of RNA metabolic process",10.029,1.9050,0.754,0.684,"secondary metabolite biosynthesis"),
c("GO:0044281","small molecule metabolic process",15.138,4.9626,0.824,0.300,"secondary metabolite biosynthesis"),
c("GO:0016053","organic acid biosynthetic process",4.171,1.9050,0.674,0.535,"secondary metabolite biosynthesis"),
c("GO:0070813","hydrogen sulfide metabolic process",0.059,2.3999,0.819,0.664,"secondary metabolite biosynthesis"),
c("GO:0019748","secondary metabolic process",0.138,4.2652,0.877,0.116,"secondary metabolite biosynthesis"),
c("GO:0006629","lipid metabolic process",3.522,1.6754,0.825,0.158,"secondary metabolite biosynthesis"),
c("GO:0006082","organic acid metabolic process",9.086,4.2652,0.732,0.575,"secondary metabolite biosynthesis"),
c("GO:0044272","sulfur compound biosynthetic process",1.235,4.7077,0.738,0.139,"secondary metabolite biosynthesis"),
c("GO:1901659","glycosyl compound biosynthetic process",1.466,2.9464,0.706,0.356,"secondary metabolite biosynthesis"),
c("GO:1901657","glycosyl compound metabolic process",2.961,2.2075,0.777,0.668,"secondary metabolite biosynthesis"),
c("GO:1901615","organic hydroxy compound metabolic process",0.831,3.2601,0.921,0.016,"organic hydroxy compound metabolism"),
c("GO:0009058","biosynthetic process",31.611,2.1545,0.946,0.027,"biosynthesis"),
c("GO:0010431","seed maturation",0.003,1.9050,0.928,0.048,"seed maturation"),
c("GO:0015772","oligosaccharide transport",0.020,2.4458,0.923,0.053,"oligosaccharide transport"),
c("GO:0042744","hydrogen peroxide catabolic process",0.093,2.0500,0.918,0.054,"hydrogen peroxide catabolism"),
c("GO:0042430","indole-containing compound metabolic process",0.285,2.3791,0.841,0.060,"indole-containing compound metabolism"),
c("GO:0006575","cellular modified amino acid metabolic process",0.785,1.6092,0.837,0.264,"indole-containing compound metabolism"),
c("GO:0042435","indole-containing compound biosynthetic process",0.194,1.9464,0.799,0.234,"indole-containing compound metabolism"),
c("GO:0018208","peptidyl-proline modification",0.406,1.8780,0.882,0.063,"peptidyl-proline modification"),
c("GO:0006415","translational termination",0.201,1.8018,0.730,0.285,"peptidyl-proline modification"),
c("GO:0051259","protein oligomerization",0.188,1.5079,0.904,0.688,"peptidyl-proline modification"),
c("GO:0001558","regulation of cell growth",0.094,1.4476,0.829,0.369,"peptidyl-proline modification"),
c("GO:0019725","cellular homeostasis",1.253,1.8018,0.800,0.072,"cellular homeostasis"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="GO_terms/ReviGO_plots/CvT_up_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the treemap command documentation for all possible parameters - there are a lot more
treemap(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "Cabbage vs TO1000 Higher Expressed Genes",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
