# ATACseq WGCNA

# libraries
library(WGCNA)
library(ggplot2)
library(gprofiler2)
library(lme4)
library(lmerTest)
library(reshape2)


# some functions
se <- function(x) {
  # standard standard error function
  sd(x)/sqrt(length(x))
}

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

###############################################################################
#### STEP 1. Data input, cleaning and pre-processing


###############################################################################
#### 1.a Loading expression data

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load dba object

# this reads in the normalized counts file from limma
data <- read.csv("input files/limma normalized DAP counts.csv", row.names = 1)

# load trait data to allow subsetting
my.design <-read.csv("input files/limma normalized DAP sample data.csv", row.names = 1)
rownames(my.design) <- gsub("X", "", rownames(my.design))


###############################################################################
#### 1.b Checking data for excessive missing values and identification of outlier microarray samples

datExpr <- data
colnames(datExpr) <- gsub("X", "", colnames(datExpr))

# Adjust file so it matches format WGCNA needs
datExpr <- as.data.frame(t(datExpr))

# run this to check if there are gene outliers
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

datExpr0 <- datExpr

# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK


# PCA
x <- prcomp(datExpr0)
PCA <- NULL
PCA$PC1 <- x$x[, 1]
PCA$PC2 <- x$x[, 2]
PCA <- as.data.frame(PCA)

ggplot(PCA, aes(PC1, PC2, label = rownames(PCA))) + 
  #geom_point(size = 5, alpha = 0.75) +
  geom_text() +
  theme_bw()


# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree <- hclust(dist(datExpr0), method = "average")
sizeGrWindow(9,6)
par(cex = 0.6, mar = c(5, 5, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
cut <- 180
# Plot a line to show the cut
abline(h = cut, col = "red")
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = cut, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = cut, col = "red")


###############################################################################
#### 1.c Loading clinical trait data
datTraits <- my.design[rownames(my.design) %in% rownames(datExpr), 
                       c("Sex", "Trt", "Side")]

# Re-cluster samples
sampleTree2 <- hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
datTraitsNumeric <- datTraits
datTraitsNumeric$Sex <- as.numeric(as.factor(datTraitsNumeric$Sex))
datTraitsNumeric$Trt <- as.numeric(as.factor(datTraitsNumeric$Trt))
datTraitsNumeric$Side <- as.numeric(as.factor(datTraitsNumeric$Side))

traitColors <- numbers2colors(datTraitsNumeric, signed = TRUE)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors, 
                    groupLabels = names(datTraitsNumeric), 
                    main = "Sample dendrogram and trait heatmap")

# save data

# check rownames of both datasets are identical
identical(rownames(datExpr), rownames(datTraits))

filename <- paste("ATACseq_", format(Sys.time(), "%Y-%m-%d"), sep = "")
save(datExpr, datTraits, file = paste("dataInput_", filename, ".RData", sep = ""))

###############################################################################
#### 2.a Automatic network construction and module detection

###############################################################################
#### 2.a.1 Choosing the soft-thresholding power: analysis of network topology

## run on biocluster as WGCNAscript2a1.R 


###############################################################################
#### 2.a.2 One-step network construction and module detection

## run on biocluster as WGCNAscript2a2.R 


###############################################################################
#### STEP 3. Relating modules to external clinical traits

### ALWAYS do this section ####

# load expression and trait file
# file <- file.choose()
file <- "input files/dataInput_ATACseq_2024-04-24.RData"
load(file)

# load results from step 2.a.2
# file <- file.choose()
file <- "input files/stickle-networkConstr-auto_dataInput_ATACseq_2024-04-24_8_30_0.15_3_0.995"
load(file)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# sample sizes
table(datTraits$Trt, datTraits$Sex)


###############################################################################
# size of modules
tableModuleColors <- table(moduleColors)
tableModuleColors 

# plot
plot(c(0, sum(tableModuleColors)), c(0, 1), type = "n", xlab = "", ylab = "", yaxt = "n")
add_x <- 0
color_order <- gsub("ME", "", names(MEs))
for (i in 1:dim(tableModuleColors)) {
  rect(add_x, 0, tableModuleColors[color_order[i]] + add_x, 1, 
       col = color_order[i])
  add_x <- add_x + tableModuleColors[color_order[i]]
}



###############################################################################
# Plot the relationships among the eigengenes

plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 4.5), marHeatmap = c(3, 4, 1, 2), 
                      cex.lab = 0.8, xLabelsAngle = 90, plotAdjacency = FALSE)


###############################################################################
# Eigengene-trait linear modeling

pedigree <- read.csv("input files/Pedigree.csv", row.names = 1, colClasses = "character")
rownames(pedigree) <- gsub("X", "", rownames(pedigree))

# this merge step was fucking shit up, fixed now
datTraits <- merge(datTraits, pedigree, by = "row.names", sort = FALSE)

# creating variables for linear modeling on MEs.
Sex <- datTraits$Sex
CareTrt <- datTraits$Trt
SIDE.RNA <- datTraits$Side
Dad <- datTraits$Dad
Family <- datTraits$Family

model <- as.formula("MEs[, i] ~  CareTrt * Sex + (1 | Dad/CareTrt/SIDE.RNA)")
num_terms <- length(attr(terms(model), "term.labels"))

significance <- matrix(NA, ncol(MEs), num_terms - 1)


for(i in 1:ncol(MEs)) {
  l <- lmer(model, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), 
                                         check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) # gets rid of singular warning 
  significance[i, ] <- anova(l, type = "II")[, 6]
  # hist(resid(l), main = paste("Residuals", colnames(MEs)[i]), xlab = "")
}
rownames(significance) <- gsub("ME", "", colnames(MEs))
colnames(significance) <- rownames(anova(l, type = "II"))

significance

write.csv(significance, "Parental_Care_LM Accessibility.csv")

# convert results to significance stars
sigStars <- matrix(NA, nrow(significance), ncol(significance))
sigStars[which(significance >= 0.1)] <- ""
sigStars[which(significance < 0.1)] <- "#"
sigStars[which(significance < 0.05)] <- "*"
sigStars[which(significance < 0.01)] <- "**"
sigStars[which(significance < 0.001)] <- "***"
rownames(sigStars) <- rownames(significance)
colnames(sigStars) <- colnames(significance)
sigStars


#########################
# eFDR

num_perm <- 10000

model <- as.formula("newMEs[, i] ~  CareTrt * Sex + (1 | Dad/CareTrt/SIDE.RNA)")
significance_eFDR <- matrix(0, ncol(MEs), num_terms - 1)
rownames(significance_eFDR) <- rownames(significance)
colnames(significance_eFDR) <- colnames(significance)

for (z in 1:num_perm) {
  newMEs <- MEs[sample(nrow(MEs)), ]
  
  print(z)
  
  significance_perm <- matrix(NA, ncol(newMEs), num_terms - 1)
  
  for(i in 1:ncol(newMEs)) {
    l <- lmer(model, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), 
                                           check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) # gets rid of singular warning 
    significance_perm[i, ] <- anova(l, type = "II")[, 6]
  }
  rownames(significance_perm) <- gsub("ME", "", colnames(newMEs))
  colnames(significance_perm) <- rownames(anova(l, type = "II"))
  
  significance_eFDR <- significance_eFDR + (significance_perm <= significance)
}

significance_eFDR <- significance_eFDR / num_perm

significance_eFDR

write.csv(significance_eFDR, "significance_eFDR Accessibility.csv")

significance_eFDR <- read.csv("significance_eFDR Accessibility.csv")


##################### figure
l <- t(MEs)
l <- l[rownames(l) != "MEgrey", ]
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")

Data <- cbind(MEs, CareTrt, Sex)
summaryTable.X <- melt(aggregate(Data[1:ncol(MEs)], by = Data[c("CareTrt", "Sex")],
                                 function(x) mean(x)))
summaryTable.SE <- melt(aggregate(Data[1:ncol(MEs)], by = Data[c("CareTrt", "Sex")],
                                  function(x) se(x)))
summaryTable <- cbind(summaryTable.X, summaryTable.SE[, 4])
names(summaryTable)[3] <- "module"
names(summaryTable)[4] <- "mean"
names(summaryTable)[5] <- "SE"

# remove grey module
summaryTable <- summaryTable[summaryTable$module != "MEgrey", ]

# make levels same as in correlation plot
summaryTable$module <- factor(summaryTable$module, levels = rev(MEs_dend$labels[MEs_dend$order]), ordered = TRUE)

pdf("Accessibility Eigengene linear models.pdf", height = 6, width = 6)
fontsize <- 18
colors <- c("#53B7E0", "#E0A473")
module <- ggplot(summaryTable, 
                 aes(x = module, y = mean, color = CareTrt, shape = Sex)) + 
  geom_point(position = position_dodge(width = 0.5), size = 5) + geom_hline(yintercept = 0) +
  coord_flip() +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE, color = CareTrt), width = 0, linewidth = 1, 
                position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("mean +/- SE Eigenpeak Accessibility") + 
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) + 
  scale_shape_discrete(labels = c("Female", "Male")) + 
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize))
module
dev.off()


# individual modules with significant effects
fontsize <- 18
GreenX <- ggplot(summaryTable[summaryTable$module == "MEgreen", ], 
                 aes(x = module, y = mean, color = CareTrt, shape = Sex)) + 
  geom_point(position = position_dodge(width = 1), size = 10) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  ggtitle("green module eigenpeak accessibility") + 
  ylim(c(ggplot_build(module)$layout$panel_scales_y[[1]]$range$range)) + 
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE, colour = CareTrt), width = 0, size = 2, position = position_dodge(width = 1)) + 
  theme_bw() +
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) + 
  scale_shape_discrete(labels = c("Female", "Male")) + 
  ylab("mean +/- SE Eigenpeak Accessibility") + 
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = fontsize))
GreenX

# individual modules with significant effects
fontsize <- 18
LightcyanX <- ggplot(summaryTable[summaryTable$module == "MElightcyan", ], 
                 aes(x = module, y = mean, color = CareTrt, shape = Sex)) + 
  geom_point(position = position_dodge(width = 1), size = 10) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  ggtitle("lightcyan module eigenpeak accessibility") + 
  ylim(c(ggplot_build(module)$layout$panel_scales_y[[1]]$range$range)) + 
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE, colour = CareTrt), width = 0, size = 2, position = position_dodge(width = 1)) + 
  theme_bw() +
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) + 
  scale_shape_discrete(labels = c("Female", "Male")) + 
  ylab("mean +/- SE Eigenpeak Accessibility") + 
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = fontsize))
LightcyanX


################################################
#### GO analyses

Peaks <- read.table("input files/ATACseq_ShiftExtend_2025_01_21_annotated_output.txt", sep = "\t", header = TRUE)
names(Peaks)[1] <- "PeakID"

set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e111_eg58_p18/")

# test
DAP_turquoise <- Peaks[Peaks$PeakID %in% colnames(datExpr[moduleColors == "turquoise"]), ]

gost.res <- gost(unique(DAP_turquoise$Entrez.ID), organism = "gaculeatus", user_threshold = 0.05, 
                 sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
print(gostplot(gost.res))
gost.res$result[-c(ncol(gost.res$result))]

dir.create(file.path("GO_Modules"), showWarnings = FALSE)

# now loop
for (color in unique(moduleColors)[unique(moduleColors) != "grey"]) {
  
  module <- Peaks[Peaks$PeakID %in% colnames(datExpr[moduleColors == color]), ]
  
  gost.res <- gost(unique(module$Entrez.ID), organism = "gaculeatus", user_threshold = 0.05, 
                   sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
  
  if (!is.null(gost.res)) {
    # gostplot(gost.res)
    write.csv(gost.res$result[-c(ncol(gost.res$result))], 
              paste0("GO_Modules/GOenrichment_gprofiler2_", color, ".csv"))
  }
  else {
    write.csv(NULL, paste0("GO_Modules/GOenrichment_gprofiler2_", color, ".csv"))
  }
}


########################################################################
# merge files produced by GO analysis

file_list <- list.files(path = "GO_Modules/", pattern = "GOenrichment_gprofiler2_")
empty <- file_list[file.size(paste0("GO_Modules/", file_list)) < 5]
file_list <- file_list[file.size(paste0("GO_Modules/", file_list)) > 5]
# file_list <- file_list[!(file_list %in% "GOenrichment_gprofiler2_grey.csv")] 

module <- gsub("GOenrichment_gprofiler2_", "", gsub(".csv", "", file_list[1]))

GOs <- read.csv(paste0("GO_Modules/", file_list[1]), row.names = 1)[c("p_value", "term_id", "term_name", "source")]
names(GOs)[1] <- module

for (i in 2:(length(file_list))) {
  file <- paste0("GO_Modules/", file_list[i])
  temp <- read.csv(file, row.names = 1)[c("p_value", "term_id", "term_name", "source")]
  module <- gsub("GOenrichment_gprofiler2_", "", gsub(".csv", "", file_list[i]))
  GOs <- merge(GOs, temp, by = c("term_id", "term_name", "source"), all = TRUE)
  names(GOs)[i + 3] <- module
}

rownames(GOs) <- paste0(GOs$term_id, ", ", GOs$term_name)

GOs <- GOs[, !names(GOs) %in% c("term_id", "term_name")]

### add empty modules back in
for (i in 1:length(empty)) {
  GOs[, paste0(gsub("GOenrichment_gprofiler2_", "", gsub(".csv", "", empty[i])))] <- NA
}

write.csv(GOs, "ATAC WGCNA Module GO.csv")

GOs_mat_BP <- as.matrix(GOs[GOs$source == "GO:BP", 2:dim(GOs)[2]])
GOs_mat_CC <- as.matrix(GOs[GOs$source == "GO:CC", 2:dim(GOs)[2]])
GOs_mat_MF <- as.matrix(GOs[GOs$source == "GO:MF", 2:dim(GOs)[2]])


########################################################################
# distribution of module genes across chromosomes
chromInfo <- read.csv("input files/seqnames.csv", row.names = 1)
chromInfo2 <- as.data.frame(chromInfo[rownames(chromInfo) %in% colnames(datExpr), ])
names(chromInfo2) <- "chromosome"
chromInfo2$module <- moduleColors
rownames(chromInfo2) <- colnames(datExpr)

chromInfo2_wo_scaffold <- chromInfo2[-(grep("scaffold", chromInfo2$chromosome)), ]

moduleChromosome <- table(chromInfo2_wo_scaffold$module, chromInfo2_wo_scaffold$chromosome)
moduleChromosome
round(moduleChromosome / rowSums(moduleChromosome) * 100, 2)
moduleChromosome <- as.data.frame(moduleChromosome)

ggplot(moduleChromosome, aes(x = Var1, y = Var2, size = Freq)) + geom_point() + theme_bw()


########################################################################
# distribution of DAPs in modules

sexDAP <- read.csv("input files/DAPs_FvsM.csv", row.names = 1)

Peaks <- read.table("input files/ATACseq_ShiftExtend_2025_01_21_annotated_output.txt", sep = "\t", header = TRUE)
names(Peaks)[1] <- "PeakID"
DAP_sex <- Peaks[Peaks$PeakID %in% rownames(sexDAP[sexDAP$adj.P.Val <= 0.05, ]), ]

# first compare distribution with modules:
DAPChromosome <- table(DAP_sex$Chr)
DAPChromosome
DAPChromosome <- as.data.frame(DAPChromosome)
DAPChromosome_wo_scaffold <- DAPChromosome[-(grep("scaffold", DAPChromosome$Var1)), ]

ggplot(DAPChromosome_wo_scaffold, aes(x = 1, y = Var1, size = Freq)) + geom_point() + theme_bw()


chromInfo2$DAP <- rownames(chromInfo2) %in% DAP_sex$PeakID 

table(chromInfo2$module, chromInfo2$DAP)

