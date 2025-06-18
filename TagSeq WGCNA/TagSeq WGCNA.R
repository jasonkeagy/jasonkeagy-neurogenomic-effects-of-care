## TagSeq WGCNA
# last run 2025 05 27
# written by Jason Keagy

# R version 4.5.0 (2025-04-11)

# libraries
library(WGCNA) # version 1.73
library(ggplot2) # version 3.5.2
library(reshape2) # version 1.4.4
library(lme4) # version 1.1-37
library(lmerTest) # version 3.1-3


# some functions
se <- function(x) {
  # standard standard error function
  sd(x)/sqrt(length(x))
}

meanplus <- function(x) { 
  # returns mean plus standard error
  mean(x) + se(x) 
}

meanminus <- function(x) { 
  # returns mean minus standard error
  mean(x) - se(x) 
}


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

###############################################################################
#### STEP 1. Data input, cleaning and pre-processing


###############################################################################
#### 1.a Loading expression data

# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# this reads in the normalized counts file from limma
data <- read.csv("input files/Limma_normalized.counts.csv", row.names = 1)

# load trait data to allow subsetting
my.design <-read.csv("input files/Limma_normalized.counts_sampledata.csv", row.names = 1)


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

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree <- hclust(dist(datExpr0), method = "average")
sizeGrWindow(9,6)
par(cex = 0.6, mar = c(5, 5, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 80, col = "red")
# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 105, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)


sampleTree <- hclust(dist(datExpr), method = "average")
sizeGrWindow(9,6)
par(cex = 0.6, mar = c(5, 5, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 80, col = "red")


###############################################################################
#### 1.c Loading clinical trait data

datTraits <- my.design[my.design$sampleID %in% gsub("X", "", rownames(datExpr)), c("Tank", "Trt", "ExpCon", "SIDE.RNA", "Sex1")]

# Re-cluster samples
sampleTree2 <- hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
datTraitsNumeric <- datTraits
datTraitsNumeric$Trt <- as.numeric(as.factor(datTraitsNumeric$Trt))
datTraitsNumeric$ExpCon <- as.numeric(as.factor(datTraitsNumeric$ExpCon))
datTraitsNumeric$SIDE.RNA <- as.numeric(as.factor(datTraitsNumeric$SIDE.RNA))
datTraitsNumeric$Sex1 <- as.numeric(as.factor(datTraitsNumeric$Sex1))

traitColors <- numbers2colors(datTraitsNumeric[c("Trt", "ExpCon", "SIDE.RNA", "Sex1")], 
  signed = TRUE)

# Plot the sample dendrogram and the colors underneath.
sizeGrWindow(9, 6)
plotDendroAndColors(sampleTree2, traitColors, 
	groupLabels = names(datTraitsNumeric[c("Trt", "ExpCon", "SIDE.RNA", "Sex1")]), 
	main = "Sample dendrogram and trait heatmap")

# save data

# check rownames of both datasets are identical
identical(rownames(datExpr), rownames(datTraits))

filename <- paste("ParentalCare_", format(Sys.time(), "%Y-%m-%d"), sep = "")
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

# load expression and trait file
file <- "input files/dataInput_ParentalCare_2019-08-27.RData"
load(file)

# load results from step 2.a.2
file <- "input files/stickle-networkConstr-auto_dataInput_ParentalCare_2019-08-27_8_30_0.15_3_0.995"
load(file)

my.design <- read.csv("Limma_normalized.counts_sampledata.csv", row.names = 1)
my.design <- my.design[rownames(my.design) %in% rownames(datTraits), ]


# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)


# sample sizes
table(datTraits$Trt, datTraits$ExpCon, datTraits$Sex1)



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

par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 4.5), marHeatmap = c(3, 4, 1, 2), 
  cex.lab = 0.8, xLabelsAngle = 90, plotAdjacency = FALSE)


###############################################################################
# Eigengene-trait linear modeling

# creating variables for linear modeling on MEs.
Sex <- datTraits$Sex1
CareTrt <- datTraits$Trt
PredExp <- datTraits$ExpCon
SIDE.RNA <- datTraits$SIDE.RNA
Dad <- as.factor(substr(datTraits$Tank, 4, 6))
Order <- as.numeric(substr(rownames(datTraits), 4, 5))

model <- as.formula("MEs[, i] ~  CareTrt * PredExp * Sex + (1 | Dad/CareTrt/PredExp/SIDE.RNA)")
num_terms <- length(attr(terms(model), "term.labels"))

significance <- matrix(NA, ncol(MEs), num_terms - 1)

# quartz(width = 9, height = 6)
# par(mfrow = c(ncol(MEs) / 4, ncol(MEs) / 3))
for(i in 1:ncol(MEs)) {
  l <- lmer(model, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5), 
    check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) # gets rid of singular warning 
  significance[i, ] <- anova(l, type = "II")[, 6]
  # hist(resid(l), main = paste("Residuals", colnames(MEs)[i]), xlab = "")
}
rownames(significance) <- gsub("ME", "", colnames(MEs))
colnames(significance) <- rownames(anova(l, type = "II"))

significance

write.csv(significance, "Parental_Care_LM.csv")

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

Sex <- datTraits$Sex1
CareTrt <- datTraits$Trt
PredExp <- datTraits$ExpCon
SIDE.RNA <- datTraits$SIDE.RNA
Dad <- as.factor(substr(datTraits$Tank, 4, 6))

model <- as.formula("newMEs[, i] ~  CareTrt * PredExp * Sex + (1 | Dad/CareTrt/PredExp/SIDE.RNA)")
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

write.csv(significance_eFDR, "significance_eFDR.csv")


##################### figure
l <- t(MEs)
l <- l[rownames(l) != "MEgrey", ]
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")
plot(MEs_dend)

Data <- cbind(MEs, CareTrt, Sex, PredExp)
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

pdf("Eigengene linear models.pdf", height = 6, width = 6)
fontsize <- 18
colors <- c("#53B7E0", "#E0A473")
module <- ggplot(summaryTable[summaryTable$module != "grey", ], 
                 aes(x = module, y = mean, color = CareTrt, shape = Sex)) + 
  geom_point(position = position_dodge(width = 0.5), size = 5) + geom_hline(yintercept = 0) +
  coord_flip() +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE, color = CareTrt), width = 0, linewidth = 1, 
      position = position_dodge(width = 0.5)) + 
  theme_bw() + 
  ylab("mean +/- SE Eigengene Expression") + 
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) + 
  scale_shape_discrete(labels = c("Female", "Male")) + 
  theme(axis.text = element_text(size = (fontsize * 0.75)),
      axis.title = element_text(size = fontsize), 
      legend.text = element_text(size = (fontsize * 0.75)),
      legend.title = element_text(size = fontsize))
module
dev.off()


# individual module with significant open field assay effect

l <- t(MEs)
l <- l[rownames(l) != "MEgrey", ]
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")
plot(MEs_dend)

Data <- cbind(MEs, CareTrt, Sex, PredExp)
summaryTable.X <- melt(aggregate(Data[1:ncol(MEs)], by = Data[c("PredExp", "Sex")],
                                 function(x) mean(x)))
summaryTable.SE <- melt(aggregate(Data[1:ncol(MEs)], by = Data[c("PredExp", "Sex")],
                                  function(x) se(x)))
summaryTable <- cbind(summaryTable.X, summaryTable.SE[, 4])
names(summaryTable)[3] <- "module"
names(summaryTable)[4] <- "mean"
names(summaryTable)[5] <- "SE"

# remove grey module
summaryTable <- summaryTable[summaryTable$module != "MEgrey", ]

# make levels same as in correlation plot
summaryTable$module <- factor(summaryTable$module, levels = rev(MEs_dend$labels[MEs_dend$order]), ordered = TRUE)

fontsize <- 18
colors <- c("purple", "gold")
GreenyellowX <- ggplot(summaryTable[summaryTable$module == "MEgreenyellow", ], 
    aes(x = module, y = mean, color = PredExp, shape = Sex)) + 
  geom_point(position = position_dodge(width = 1), size = 10) + 
  geom_hline(yintercept = 0) +
  coord_flip() +
  ggtitle("greenyellow module eigengene expression") + 
  ylim(c(ggplot_build(module)$layout$panel_scales_y[[1]]$range$range)) + 
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE, colour = PredExp), width = 0, linewidth = 2, position = position_dodge(width = 1)) + 
  theme_bw() +
  scale_color_discrete(type = colors, name = "Open Assay\nExperiment", labels = c("Baseline Control", "Experimental")) + 
  scale_shape_discrete(labels = c("Female", "Male")) + 
  ylab("mean +/- SE Eigengene Expression") + 
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = fontsize))
GreenyellowX

