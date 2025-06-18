## gprofiler Gene Expression
# last run 2025 05 27
# written by Jason Keagy

# R version 4.5.0 (2025-04-11)

# load all the libraries

library(gprofiler2) # version 0.2.3
#library(org.Dr.eg.db) # version 3.21.0
#library(GOSemSim) # version 2.34.0
#library(ComplexHeatmap) # version 2.25.1
#library(circlize) # version 0.4.16
#library(MASS) # version 7.3-65
library(tidyverse) # version 2.0.0
library(WGCNA) # version 1.73
library(RColorBrewer) # version 1.1-3


# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


################################################################################

Day0vs10_DEGs_file <- "input files/DEGs_Day0vs10.csv"
Day0vs10_DEGs <- read.csv(Day0vs10_DEGs_file, row.names = 1)
Day0vs10_DEGs <- Day0vs10_DEGs[Day0vs10_DEGs$adj.P.Val < 0.05, ]

# write.table(rownames(Day0vs10_DEGs), "DEGs_Day0vs10.txt", row.names = FALSE, quote = FALSE, col.names = FALSE)

set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e111_eg58_p18/")

gost.res <- gost(rownames(Day0vs10_DEGs), organism = "gaculeatus", user_threshold = 0.05,
                 sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
if (!is.null(gost.res)) {
  print(gostplot(gost.res))
  write.csv(gost.res$result[-c(ncol(gost.res$result))], "GO_Day0vs10.csv")
} else {
  write.csv(NULL, "GO_Day0vs10.csv")
}


FvsM_DEGs_file <- "input files/DEGs_FvsM.csv"
FvsM_DEGs <- read.csv(FvsM_DEGs_file, row.names = 1)
FvsM_DEGs <- FvsM_DEGs[FvsM_DEGs$adj.P.Val < 0.05, ]

gost.res <- gost(rownames(FvsM_DEGs), organism = "gaculeatus", user_threshold = 0.05, 
                 sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
if (!is.null(gost.res)) {
  print(gostplot(gost.res))
  write.csv(gost.res$result[-c(ncol(gost.res$result))], "GO_FvsM.csv")
} else {
  write.csv(NULL, "GO_FvsM.csv")
}


CvsE_DEGs_file <- "input files/DEGs_CvsE.csv"
CvsE_DEGs <- read.csv(CvsE_DEGs_file, row.names = 1)
CvsE_DEGs <- CvsE_DEGs[CvsE_DEGs$adj.P.Val < 0.05, ]

gost.res <- gost(rownames(CvsE_DEGs), organism = "gaculeatus", user_threshold = 0.05, 
                 sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
if (!is.null(gost.res)) {
  print(gostplot(gost.res))
  write.csv(gost.res$result[-c(ncol(gost.res$result))], "GO_CvsE.csv")
} else {
  write.csv(NULL, "GO_CvsE.csv")
}


GO_Day0vs10 <- read.csv("GO_Day0vs10.csv", row.names = 1)

GO_Day0vs10$fullname <- paste(GO_Day0vs10$term_id, GO_Day0vs10$term_name)


ggplot(GO_Day0vs10, aes(x = recall, y = fullname, fill = -log10(p_value))) + 
  geom_point(size = 4, pch = 21, color = "black") + facet_grid (source ~ ., scales = "free", space = "free") +
  scale_fill_gradientn(name = "-log10(p)", colors = brewer.pal(5, "BuPu"), limits = c(1,10)) +
  xlab("Gene ratios") + ylab("") +
  xlim(c(0, 0.6)) +
  theme_bw()


################################################################################

# WGCNA modules

# load expression and trait file
file1 <- "input files/dataInput_ParentalCare_2019-08-27.RData"
load(file1)
file2 <- "input files/stickle-networkConstr-auto_dataInput_ParentalCare_2019-08-27_8_30_0.15_3_0.995"
load(file2)


# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# dendrogram of module eigengenes
l <- t(MEs)
l <- l[rownames(l) != "MEgrey", ]
rownames(l) <- gsub("ME", "", rownames(l))
MEs_dend <- hclust(as.dist(1 - cor(t(l))), method = "average")


########################################################################
# GO analysis of modules

dir.create(file.path("GO_Modules"), showWarnings = FALSE)
# repeat for each module:
for (color in gsub("ME", "", colnames(MEs))) {
  
    
  # color <- "cyan"
  
  module <- names(datExpr[moduleColors == color])
  
  gost.res <- gost(module, organism = "gaculeatus", user_threshold = 0.05, 
                   sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
  
  if (!is.null(gost.res)) {
    gostplot(gost.res)
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

GOs <- GOs[, c("source", rownames(l))]

write.csv(GOs, "Combined WGCNA GO terms.csv")




