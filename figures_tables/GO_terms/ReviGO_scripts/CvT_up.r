

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
revigo.data <- rbind(c("GO:0006790","sulfur compound metabolic process",1.822,6.2111,0.900,0.000,"sulfur compound metabolism"),
c("GO:0009625","response to insect",0.001,3.9075,0.837,0.000,"response to insect"),
c("GO:0009628","response to abiotic stimulus",0.571,2.0603,0.804,0.448,"response to insect"),
c("GO:0001101","response to acid chemical",0.124,2.4724,0.788,0.563,"response to insect"),
c("GO:0031668","cellular response to extracellular stimulus",0.436,1.5889,0.742,0.575,"response to insect"),
c("GO:1901700","response to oxygen-containing compound",0.503,2.4794,0.770,0.298,"response to insect"),
c("GO:0010438","cellular response to sulfur starvation",0.000,2.7077,0.800,0.337,"response to insect"),
c("GO:0009611","response to wounding",0.127,3.7593,0.798,0.427,"response to insect"),
c("GO:0006979","response to oxidative stress",0.575,1.5752,0.779,0.498,"response to insect"),
c("GO:0006970","response to osmotic stress",0.082,3.7593,0.781,0.195,"response to insect"),
c("GO:0006282","regulation of DNA repair",0.057,1.4633,0.674,0.404,"response to insect"),
c("GO:0010200","response to chitin",0.004,1.3566,0.807,0.441,"response to insect"),
c("GO:0009414","response to water deprivation",0.022,3.4081,0.758,0.662,"response to insect"),
c("GO:0006950","response to stress",4.575,2.1281,0.775,0.441,"response to insect"),
c("GO:0007165","signal transduction",6.621,1.3337,0.602,0.637,"response to insect"),
c("GO:0010035","response to inorganic substance",0.317,1.9629,0.776,0.611,"response to insect"),
c("GO:0050896","response to stimulus",12.210,6.4895,0.970,0.000,"response to stimulus"),
c("GO:1901615","organic hydroxy compound metabolic process",0.831,3.7593,0.919,0.016,"organic hydroxy compound metabolism"),
c("GO:0009058","biosynthetic process",31.611,1.5964,0.945,0.027,"biosynthesis"),
c("GO:0051259","protein oligomerization",0.188,2.2410,0.925,0.032,"protein oligomerization"),
c("GO:0001558","regulation of cell growth",0.094,1.5350,0.813,0.368,"protein oligomerization"),
c("GO:0019758","glycosinolate biosynthetic process",0.001,4.7620,0.654,0.038,"glycosinolate biosynthesis"),
c("GO:1901576","organic substance biosynthetic process",30.365,2.1165,0.788,0.462,"glycosinolate biosynthesis"),
c("GO:0032774","RNA biosynthetic process",10.925,3.0802,0.706,0.349,"glycosinolate biosynthesis"),
c("GO:0006575","cellular modified amino acid metabolic process",0.785,1.5257,0.829,0.264,"glycosinolate biosynthesis"),
c("GO:0044283","small molecule biosynthetic process",5.677,3.0282,0.687,0.655,"glycosinolate biosynthesis"),
c("GO:0051252","regulation of RNA metabolic process",10.029,2.3789,0.697,0.684,"glycosinolate biosynthesis"),
c("GO:0044281","small molecule metabolic process",15.138,5.2441,0.814,0.117,"glycosinolate biosynthesis"),
c("GO:0030003","cellular cation homeostasis",0.282,1.3196,0.813,0.107,"glycosinolate biosynthesis"),
c("GO:0016053","organic acid biosynthetic process",4.171,3.1139,0.644,0.535,"glycosinolate biosynthesis"),
c("GO:0042430","indole-containing compound metabolic process",0.285,2.6633,0.829,0.162,"glycosinolate biosynthesis"),
c("GO:0042435","indole-containing compound biosynthetic process",0.194,1.7339,0.784,0.234,"glycosinolate biosynthesis"),
c("GO:0055114","oxidation-reduction process",15.060,1.8417,0.814,0.415,"glycosinolate biosynthesis"),
c("GO:0006082","organic acid metabolic process",9.086,5.2441,0.715,0.575,"glycosinolate biosynthesis"),
c("GO:1901659","glycosyl compound biosynthetic process",1.466,3.8463,0.688,0.433,"glycosinolate biosynthesis"),
c("GO:0044272","sulfur compound biosynthetic process",1.235,3.4551,0.751,0.493,"glycosinolate biosynthesis"),
c("GO:1901657","glycosyl compound metabolic process",2.961,3.1193,0.770,0.668,"glycosinolate biosynthesis"),
c("GO:0015772","oligosaccharide transport",0.020,1.7671,0.923,0.041,"oligosaccharide transport"),
c("GO:0007154","cell communication",7.219,2.2410,0.924,0.049,"cell communication"),
c("GO:0007568","aging",0.088,1.6958,0.893,0.052,"aging"),
c("GO:0010431","seed maturation",0.003,1.7498,0.903,0.430,"aging"),
c("GO:0090693","plant organ senescence",0.006,1.7861,0.901,0.418,"aging"),
c("GO:0018208","peptidyl-proline modification",0.406,1.8417,0.884,0.063,"peptidyl-proline modification"),
c("GO:1901617","organic hydroxy compound biosynthetic process",0.383,1.9931,0.859,0.097,"organic hydroxy compound biosynthesis"));revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006790","sulfur compound metabolic process",1.822,6.2111,0.900,0.000,"sulfur compound metabolism"),
c("GO:0009625","response to insect",0.001,3.9075,0.837,0.000,"response to insect"),
c("GO:0009628","response to abiotic stimulus",0.571,2.0603,0.804,0.448,"response to insect"),
c("GO:0001101","response to acid chemical",0.124,2.4724,0.788,0.563,"response to insect"),
c("GO:0031668","cellular response to extracellular stimulus",0.436,1.5889,0.742,0.575,"response to insect"),
c("GO:1901700","response to oxygen-containing compound",0.503,2.4794,0.770,0.298,"response to insect"),
c("GO:0010438","cellular response to sulfur starvation",0.000,2.7077,0.800,0.337,"response to insect"),
c("GO:0009611","response to wounding",0.127,3.7593,0.798,0.427,"response to insect"),
c("GO:0006979","response to oxidative stress",0.575,1.5752,0.779,0.498,"response to insect"),
c("GO:0006970","response to osmotic stress",0.082,3.7593,0.781,0.195,"response to insect"),
c("GO:0006282","regulation of DNA repair",0.057,1.4633,0.674,0.404,"response to insect"),
c("GO:0010200","response to chitin",0.004,1.3566,0.807,0.441,"response to insect"),
c("GO:0009414","response to water deprivation",0.022,3.4081,0.758,0.662,"response to insect"),
c("GO:0006950","response to stress",4.575,2.1281,0.775,0.441,"response to insect"),
c("GO:0007165","signal transduction",6.621,1.3337,0.602,0.637,"response to insect"),
c("GO:0010035","response to inorganic substance",0.317,1.9629,0.776,0.611,"response to insect"),
c("GO:0050896","response to stimulus",12.210,6.4895,0.970,0.000,"response to stimulus"),
c("GO:1901615","organic hydroxy compound metabolic process",0.831,3.7593,0.919,0.016,"organic hydroxy compound metabolism"),
c("GO:0009058","biosynthetic process",31.611,1.5964,0.945,0.027,"biosynthesis"),
c("GO:0051259","protein oligomerization",0.188,2.2410,0.925,0.032,"protein oligomerization"),
c("GO:0001558","regulation of cell growth",0.094,1.5350,0.813,0.368,"protein oligomerization"),
c("GO:0019758","glycosinolate biosynthetic process",0.001,4.7620,0.654,0.038,"glycosinolate biosynthesis"),
c("GO:1901576","organic substance biosynthetic process",30.365,2.1165,0.788,0.462,"glycosinolate biosynthesis"),
c("GO:0032774","RNA biosynthetic process",10.925,3.0802,0.706,0.349,"glycosinolate biosynthesis"),
c("GO:0006575","cellular modified amino acid metabolic process",0.785,1.5257,0.829,0.264,"glycosinolate biosynthesis"),
c("GO:0044283","small molecule biosynthetic process",5.677,3.0282,0.687,0.655,"glycosinolate biosynthesis"),
c("GO:0051252","regulation of RNA metabolic process",10.029,2.3789,0.697,0.684,"glycosinolate biosynthesis"),
c("GO:0044281","small molecule metabolic process",15.138,5.2441,0.814,0.117,"glycosinolate biosynthesis"),
c("GO:0030003","cellular cation homeostasis",0.282,1.3196,0.813,0.107,"glycosinolate biosynthesis"),
c("GO:0016053","organic acid biosynthetic process",4.171,3.1139,0.644,0.535,"glycosinolate biosynthesis"),
c("GO:0042430","indole-containing compound metabolic process",0.285,2.6633,0.829,0.162,"glycosinolate biosynthesis"),
c("GO:0042435","indole-containing compound biosynthetic process",0.194,1.7339,0.784,0.234,"glycosinolate biosynthesis"),
c("GO:0055114","oxidation-reduction process",15.060,1.8417,0.814,0.415,"glycosinolate biosynthesis"),
c("GO:0006082","organic acid metabolic process",9.086,5.2441,0.715,0.575,"glycosinolate biosynthesis"),
c("GO:1901659","glycosyl compound biosynthetic process",1.466,3.8463,0.688,0.433,"glycosinolate biosynthesis"),
c("GO:0044272","sulfur compound biosynthetic process",1.235,3.4551,0.751,0.493,"glycosinolate biosynthesis"),
c("GO:1901657","glycosyl compound metabolic process",2.961,3.1193,0.770,0.668,"glycosinolate biosynthesis"),
c("GO:0015772","oligosaccharide transport",0.020,1.7671,0.923,0.041,"oligosaccharide transport"),
c("GO:0007154","cell communication",7.219,2.2410,0.924,0.049,"cell communication"),
c("GO:0007568","aging",0.088,1.6958,0.893,0.052,"aging"),
c("GO:0010431","seed maturation",0.003,1.7498,0.903,0.430,"aging"),
c("GO:0090693","plant organ senescence",0.006,1.7861,0.901,0.418,"aging"),
c("GO:0018208","peptidyl-proline modification",0.406,1.8417,0.884,0.063,"peptidyl-proline modification"),
c("GO:1901617","organic hydroxy compound biosynthetic process",0.383,1.9931,0.859,0.097,"organic hydroxy compound biosynthesis"));

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
