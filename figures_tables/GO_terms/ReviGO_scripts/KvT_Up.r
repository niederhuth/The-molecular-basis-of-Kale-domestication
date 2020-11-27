

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
revigo.data <- rbind(c("GO:0010035","response to inorganic substance",0.317,4.5031,0.826,0.000,"response to inorganic substance"),
c("GO:0009625","response to insect",0.001,2.8060,0.866,0.478,"response to inorganic substance"),
c("GO:0001101","response to acid chemical",0.124,2.1385,0.834,0.542,"response to inorganic substance"),
c("GO:0031668","cellular response to extracellular stimulus",0.436,1.6044,0.803,0.491,"response to inorganic substance"),
c("GO:0009611","response to wounding",0.127,2.6760,0.843,0.297,"response to inorganic substance"),
c("GO:0006970","response to osmotic stress",0.082,2.1969,0.826,0.669,"response to inorganic substance"),
c("GO:0009415","response to water",0.026,2.2850,0.804,0.482,"response to inorganic substance"),
c("GO:0006950","response to stress",4.575,2.5008,0.837,0.418,"response to inorganic substance"),
c("GO:0010363","regulation of plant-type hypersensitive response",0.001,1.7091,0.775,0.318,"response to inorganic substance"),
c("GO:0080027","response to herbivore",0.000,3.7283,0.869,0.203,"response to inorganic substance"),
c("GO:0015979","photosynthesis",0.183,9.6003,0.931,0.000,"photosynthesis"),
c("GO:0044419","interspecies interaction between organisms",0.262,3.4858,0.969,0.000,"interspecies interaction between organisms"),
c("GO:0050896","response to stimulus",12.210,6.2549,0.972,0.000,"response to stimulus"),
c("GO:0006629","lipid metabolic process",3.522,6.3990,0.833,0.014,"lipid metabolism"),
c("GO:1901570","fatty acid derivative biosynthetic process",0.009,2.2742,0.834,0.243,"lipid metabolism"),
c("GO:1901568","fatty acid derivative metabolic process",0.017,5.1451,0.884,0.131,"lipid metabolism"),
c("GO:0044283","small molecule biosynthetic process",5.677,2.4873,0.761,0.655,"lipid metabolism"),
c("GO:0044281","small molecule metabolic process",15.138,9.2111,0.824,0.300,"lipid metabolism"),
c("GO:0044550","secondary metabolite biosynthetic process",0.101,2.6468,0.786,0.153,"lipid metabolism"),
c("GO:0008610","lipid biosynthetic process",2.123,2.4873,0.766,0.563,"lipid metabolism"),
c("GO:0016053","organic acid biosynthetic process",4.171,2.5929,0.721,0.393,"lipid metabolism"),
c("GO:0019748","secondary metabolic process",0.138,1.6422,0.878,0.158,"lipid metabolism"),
c("GO:0055114","oxidation-reduction process",15.060,2.1385,0.824,0.415,"lipid metabolism"),
c("GO:0006082","organic acid metabolic process",9.086,8.8447,0.751,0.575,"lipid metabolism"),
c("GO:0015977","carbon fixation",0.036,2.6422,0.880,0.139,"lipid metabolism"),
c("GO:0009657","plastid organization",0.024,1.6977,0.939,0.022,"plastid organization"),
c("GO:0031032","actomyosin structure organization",0.059,1.3212,0.869,0.407,"plastid organization"),
c("GO:0050898","nitrile metabolic process",0.000,2.6220,0.930,0.027,"nitrile metabolism"),
c("GO:0080028","nitrile biosynthetic process",0.000,2.7769,0.908,0.103,"nitrile metabolism"),
c("GO:0009765","photosynthesis, light harvesting",0.019,1.6870,0.937,0.040,"photosynthesis, light harvesting"),
c("GO:0042537","benzene-containing compound metabolic process",0.162,1.7091,0.911,0.047,"benzene-containing compound metabolism"),
c("GO:1901617","organic hydroxy compound biosynthetic process",0.383,1.7805,0.887,0.050,"organic hydroxy compound biosynthesis"),
c("GO:1901136","carbohydrate derivative catabolic process",0.423,2.6740,0.905,0.051,"carbohydrate derivative catabolism"),
c("GO:0006026","aminoglycan catabolic process",0.166,1.3212,0.887,0.495,"carbohydrate derivative catabolism"),
c("GO:1901615","organic hydroxy compound metabolic process",0.831,6.1701,0.935,0.055,"organic hydroxy compound metabolism"),
c("GO:0018298","protein-chromophore linkage",0.095,4.5031,0.918,0.056,"protein-chromophore linkage"),
c("GO:0006790","sulfur compound metabolic process",1.822,2.5929,0.922,0.058,"sulfur compound metabolism"),
c("GO:0010119","regulation of stomatal movement",0.004,3.5972,0.855,0.060,"regulation of stomatal movement"),
c("GO:1903792","negative regulation of anion transport",0.004,2.1385,0.820,0.174,"regulation of stomatal movement"),
c("GO:0019725","cellular homeostasis",1.253,2.8060,0.782,0.173,"regulation of stomatal movement"),
c("GO:0006355","regulation of transcription, DNA-templated",9.917,1.8908,0.755,0.377,"regulation of stomatal movement"),
c("GO:0090693","plant organ senescence",0.006,2.6468,0.880,0.062,"plant organ senescence"),
c("GO:0019098","reproductive behavior",0.019,1.8400,0.943,0.398,"plant organ senescence"),
c("GO:0007568","aging",0.088,2.2230,0.895,0.418,"plant organ senescence"),
c("GO:0010118","stomatal movement",0.006,2.5929,0.893,0.094,"stomatal movement"));

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
