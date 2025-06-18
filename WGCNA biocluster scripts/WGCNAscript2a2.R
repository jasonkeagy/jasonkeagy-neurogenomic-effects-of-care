# script for running WGCNA on biocluster

args <- commandArgs(trailingOnly = TRUE)

library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


###############################################################################
#### STEP 2. Automatic construction of the gene network and identification of modules

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# Make sure to set this correctly.
enableWGCNAThreads(nThreads = 12)

# Load the data saved in the first part
lnames <- load(file = args[1])
print(paste("file = ", args[1], sep = ""))
# The variable lnames contains the names of loaded variables
lnames


###############################################################################
#### 2.a Automatic network construction and module detection


###############################################################################
#### 2.a.1 Choosing the soft-thresholding power: analysis of network topology

powernum <- as.numeric(args[2])
print(paste("powernum = ", powernum, sep = ""))


###############################################################################
#### 2.a.2 One-step network construction and module detection

# Constructing the gene network and identifying modules is now a simple function call:
# set this up to be one of the parameters that is passed to the script.
max <- 50000 
# min_module_size <- 30
min_module_size <- as.numeric(args[3])
print(paste("min module size = ", min_module_size, sep = ""))
# merge_cut_height <- 0.25, default is actually 0.15
merge_cut_height <- as.numeric(args[4])
print(paste("merge cut height = ", merge_cut_height, sep = ""))
# deep_split <- 2 # varies between 0-4
deep_split <- as.numeric(args[5])
print(paste("deep split = ", deep_split, sep = ""))
# detectCutHeight <- 0.995
detect_cut_height <- as.numeric(args[6])
print(paste("detectCutHeight = ", detect_cut_height, sep = ""))


net <- blockwiseModules(datExpr, 
	                    maxBlockSize = max, 
	                    power = powernum, 
	                    networkType = "signed",
                        corType = "bicor", 
                        maxPOutliers = 0.05,
                        minModuleSize = min_module_size, 
	                    mergeCutHeight = merge_cut_height, 
	                    numericLabels = TRUE, 
	                    saveTOMs = TRUE,
	                    saveTOMFileBase = paste("stickleTOM_", gsub(".RData", "", args[1]), 
	                    	args[2], args[3], args[4], args[5], args[6], sep = "_"), 
	                    verbose = 3, 
	                    deepSplit = deep_split, 
	                    detectCutHeight = detect_cut_height)


# "recutBlockwiseTrees" can be used to make changes without having to redo stuff (e.g., tree cut, module membership, and module merging criteria)


# We now save the module assignment and module eigengene information necessary for subsequent analysis.
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = paste("stickle-networkConstr-auto", 
	gsub(".RData", "", args[1]), args[2], args[3], args[4], args[5], args[6], sep = "_"))


###############################################################################

# open a graphics window
pdf(paste("stickle-networkConstr-auto", gsub(".RData", "", args[1]), args[2], 
	args[3], args[4], args[5], args[6], ".pdf", sep = "_"), 6, 6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], "Module colors", 
	main = "Gene dendrogram and module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
	guideHang = 0.05)
dev.off()


