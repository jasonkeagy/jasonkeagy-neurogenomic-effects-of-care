## ATACseq TagSeq WGCNA overlap
# last run 2025 05 27
# written by Jason Keagy

# R version 4.5.0 (2025-04-11)

# libraries
library(ComplexHeatmap) # version 2.25.1
library(WGCNA) # version 1.73
library(circlize) # version 0.4.16


# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# ATACseq modules
file <- "input files/stickle-networkConstr-auto_dataInput_ATACseq_2024-04-24_8_30_0.15_3_0.995"
load(file)

# ATAC data
file <- "input files/dataInput_ATACseq_2024-04-24.RData"
load(file)

ATAC_datExpr <- datExpr

ATAC_MEs <- MEs
ATAC_moduleColors <- moduleColors
ATAC_moduleLabels <- moduleLabels


# expression modules
file <- "input files/stickle-networkConstr-auto_dataInput_ParentalCare_2019-08-27_8_30_0.15_3_0.995"
load(file)

# expression data
file <- "input files/dataInput_ParentalCare_2019-08-27.RData"
load(file)

WGCNA_datExpr <- datExpr

WGCNA_MEs <- MEs
WGCNA_moduleColors <- moduleColors
WGCNA_moduleLabels <- moduleLabels

WGCNA_moduleLabels <- data.frame(WGCNA_moduleLabels)
WGCNA_moduleLabels$gene <- rownames(WGCNA_moduleLabels)


Peaks <- read.table("input files/ATACseq_ShiftExtend_2025_01_21_annotated_output.txt", sep = "\t", header = TRUE)
names(Peaks)[1] <- "PeakID"

ATAC_moduleLabels_annot <- merge(ATAC_moduleLabels, Peaks[c("PeakID", "Entrez.ID")], by.x = "row.names", by.y = "PeakID")

unique_WGCNA_modules <- unique(WGCNA_moduleLabels$WGCNA_moduleLabels)
unique_ATAC_modules <- unique(ATAC_moduleLabels_annot$x)


overlap <- as.data.frame(matrix(NA, length(unique_WGCNA_modules), length(unique_ATAC_modules)))
rownames(overlap) <- unique(WGCNA_moduleColors)
colnames(overlap) <- unique(ATAC_moduleColors)

for (i in unique_WGCNA_modules) {
  for (j in unique_ATAC_modules) {
    #print(paste("i=", i, "j=", j))
    #print(sum(WGCNA_moduleLabels$gene[WGCNA_moduleLabels$WGCNA_moduleLabels == i] %in% 
    #            ATAC_moduleLabels_annot$Entrez.ID[ATAC_moduleLabels_annot$x == j]))
  overlap[i+1, j+1] <- sum(WGCNA_moduleLabels$gene[WGCNA_moduleLabels$WGCNA_moduleLabels == i] %in% 
                         ATAC_moduleLabels_annot$Entrez.ID[ATAC_moduleLabels_annot$x == j])
  }
}

overlap

apply(overlap, 1, FUN = function(x) {x / sum(x) * 100})
apply(overlap, 2, FUN = function(x) {x / sum(x) * 100})



overlap_p <- as.data.frame(matrix(0, length(unique_WGCNA_modules), length(unique_ATAC_modules)))
rownames(overlap_p) <- unique(WGCNA_moduleColors)
colnames(overlap_p) <- unique(ATAC_moduleColors)

num_perm <- 10000

for (k in 1:num_perm) {
  print(k)
  ATAC_moduleLabels_annot$x <- sample(ATAC_moduleLabels_annot$x)
  WGCNA_moduleLabels$WGCNA_moduleLabels <- sample(WGCNA_moduleLabels$WGCNA_moduleLabels)
  for (i in unique_WGCNA_modules) {
    for (j in unique_ATAC_modules) {
      if(overlap[i+1, j+1] >= sum(WGCNA_moduleLabels$gene[WGCNA_moduleLabels$WGCNA_moduleLabels == i] %in% 
                                 ATAC_moduleLabels_annot$Entrez.ID[ATAC_moduleLabels_annot$x == j])) {
        overlap_p[i+1, j+1] <- overlap_p[i+1, j+1] + 1
      }
    }
  }
}


overlap_pvalue <- overlap_p/num_perm

write.csv(overlap_pvalue, "Expression accessbility WGCNA overlap.csv")

overlap_pvalue <- read.csv("Expression accessbility WGCNA overlap.csv", row.names = 1)

overlap_pvalue <- overlap_pvalue[rownames(overlap_pvalue) != "grey", colnames(overlap_pvalue) != "grey"]


MEs0 = moduleEigengenes(ATAC_datExpr, ATAC_moduleColors)$eigengenes
MEs = orderMEs(MEs0)
l <- t(MEs)
rownames(l) <- gsub("ME", "", rownames(l))
l <- l[rownames(l) %in% colnames(overlap_pvalue), ]
ATAC_MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")

MEs0 = moduleEigengenes(WGCNA_datExpr, WGCNA_moduleColors)$eigengenes
MEs = orderMEs(MEs0)
l <- t(MEs)
rownames(l) <- gsub("ME", "", rownames(l))
l <- l[rownames(l) %in% rownames(overlap_pvalue), ]
WGCNA_MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")


overlap_pvalue <-  overlap_pvalue[WGCNA_MEs_dend$labels, ATAC_MEs_dend$labels]

col_fun = colorRamp2(c(0, 0.055), c("red", "white"))

Heatmap(as.matrix(overlap_pvalue), col = col_fun,
        cluster_rows = WGCNA_MEs_dend, cluster_columns = ATAC_MEs_dend,
        name = "p-value"
        )

