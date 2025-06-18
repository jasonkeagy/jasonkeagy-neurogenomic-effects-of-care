# script for running WGCNA on biocluster

args <- commandArgs(trailingOnly = TRUE)


library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


###############################################################################
#### STEP 2. Automatic construction of the gene network and identification of modules

# Load the data saved in the first part
lnames <- load(file = args[1])
print(paste("file = ", args[1], sep = ""))
# The variable lnames contains the names of loaded variables
lnames


###############################################################################
#### 2.a Automatic network construction and module detection

###############################################################################
#### 2.a.1 Choosing the soft-thresholding power: analysis of network topology

# unsigned
# Choose a set of soft-thresholding powers
powers <- c(1:20)
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "unsigned")
# Plot the results:
pdf(paste("step_2a", gsub(".RData", "", args[1]), "unsigned.pdf", sep = "_"), 10, 6)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[,3]) * sft$fitIndices[, 2], 
  xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
  main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
    labels = powers, cex = cex1, col = "red");
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()


#signed
# Choose a set of soft-thresholding powers
powers <- c(1:20)
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
pdf(paste("step_2a", gsub(".RData", "", args[1]), "signed.pdf", sep = "_"), 10, 6)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[,3]) * sft$fitIndices[, 2], 
  xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
  main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2],
    labels = powers, cex = cex1, col = "red");
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.90, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

