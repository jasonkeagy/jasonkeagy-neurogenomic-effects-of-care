## ATACseq Analysis
# last run 2025 05 27
# written by Jason Keagy

# R version 4.5.0 (2025-04-11)

# libraries
library(DiffBind) # version 3.18.0
library(limma) # version 3.64.1
library(ChIPseeker) # version 1.44.0
library(edgeR) # version 4.6.2
library(Vennerable) # version 3.1.0.900
library(gprofiler2) # version 0.2.3
library(ATACseqQC) # version 1.32.0
library(ggplot2) # version 3.5.2


# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# load dba object

ss <- dba.load('dba_ATACseq_ShiftExtend_2025_01_21', pre = 'input files/')

# Convert DBA object to SummarizedExperiment
ss_sum <- dba(ss, bSummarizedExperiment = TRUE)
ss_sum

# convert factor names to more interpretable names
names(ss_sum@colData)[c(1:2)] <- c("Sex", "Trt")
ss_sum@colData$Side <- ss$samples$Side

write.csv(seqnames(ss_sum@rowRanges), "seqnames.csv")


# ---------------------------------------------------------------------------------------------

ss_txt <- cbind(rownames(ss$peaks[[1]]), ss$peaks[[1]][c(1:3)], "+")
names(ss_txt) <- c("#PeakID", "Chr", "Start", "End", "Strand")

ss_bed <- cbind(ss$peaks[[1]][c(1:3)], rownames(ss$peaks[[1]]), NA, "+")
names(ss_bed) <- c("#Chr", "Start", "End", "PeakID", "", "Strand")

write.table(ss_txt, "ATACseq_ShiftExtend_2025_01_21.txt", sep = "\t", row.names = FALSE, quote = FALSE, eol = "\n")
write.table(ss_bed, "ATACseq_ShiftExtend_2025_01_21.bed", sep = "\t", row.names = FALSE, quote = FALSE, eol = "\n")

# annotatePeaks.pl ATACseq_ShiftExtend_2025_01_21.txt Gasterosteus_aculeatus.BROADS1.99.dna_rm.toplevel.fa -gtf Gasterosteus_aculeatus.BROADS1.99.gtf

# ---------------------------------------------------------------------------------------------


# limma

ss_sum$Group <- factor(paste0(ss_sum$Trt, "_", ss_sum$Sex))
mm <- model.matrix(~ 0 + ss_sum$Group + ss_sum$Side) 
colnames(mm) <- gsub("Group", "", gsub("ss_sum[$]", "", colnames(mm)))

d <- assay(ss_sum, "Reads")
d <- DGEList(d)

keep <- filterByExpr(d, mm, min.count = 10) # filter using model matrix specifications of cell group sizes
d <- d[keep, , keep.lib.sizes = FALSE]

d <- calcNormFactors(d)
v <- voom(d, mm, plot = TRUE)
write.csv(v$E, "Limma normalized DAP counts.csv")
write.csv(colData(ss_sum), "Limma normalized DAP sample data.csv")

fit <- lmFit(v, mm)


contrasts <- makeContrasts(
  contrast_0vs10 = (Day0_F + Day0_M) / 2 - (Day10_F + Day10_M) / 2, 
  contrast_FvsM = (Day0_F + Day10_F) / 2 - (Day0_M + Day10_M) / 2, 
  contrast_0vs10_SexInt = (Day0_M - Day10_M) - (Day0_F - Day10_F),
  levels = mm)


fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)


top.table_0vs10 <- topTable(fit2, coef = "contrast_0vs10", sort.by = "P", n = Inf)
write.csv(top.table_0vs10, "DAPs_Day0vs10.csv")

top.table_FvsM <- topTable(fit2, coef = "contrast_FvsM", sort.by = "P", n = Inf)
write.csv(top.table_FvsM, "DAPs_FvsM.csv")

top.table_0vs10_SexInt <- topTable(fit2, coef = "contrast_0vs10_SexInt", sort.by = "P", n = Inf)
write.csv(top.table_0vs10_SexInt, "DAPs_Day0vs10_SexInt.csv")


nrow(top.table_0vs10[top.table_0vs10$adj.P.Val <= 0.05, ])
nrow(top.table_FvsM[top.table_FvsM$adj.P.Val <= 0.05, ])
nrow(top.table_0vs10_SexInt[top.table_0vs10_SexInt$adj.P.Val <= 0.05, ])

nrow(top.table_0vs10[top.table_0vs10$adj.P.Val <= 0.1, ])
nrow(top.table_FvsM[top.table_FvsM$adj.P.Val <= 0.1, ])
nrow(top.table_0vs10_SexInt[top.table_0vs10_SexInt$adj.P.Val <= 0.1, ])


size = 18
vol_sexDAP <- ggplot(top.table_FvsM, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(data = top.table_FvsM[top.table_FvsM$adj.P.Val > 0.05, ],
             color = "grey", shape = 19, size = 2, alpha = 0.5) +
  geom_point(data = top.table_FvsM[top.table_FvsM$adj.P.Val <= 0.05 &
                                     top.table_FvsM$logFC > 0, ],
             color = "blue", shape = 19, size = 2, alpha = 0.5) +
  geom_point(data = top.table_FvsM[top.table_FvsM$adj.P.Val <= 0.05 &
                                     top.table_FvsM$logFC < 0, ],
             color = "red", shape = 19, size = 2, alpha = 0.5) +
  xlab(expression(log[2]~fold~change)) + 
  ylab(expression(-log[10]~(FDR))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  theme(axis.title = element_text(size = size),
        axis.text = element_text(size = size - 2),
        legend.title = element_text(size = size), 
        legend.text = element_text(size = size - 2))
vol_sexDAP

pdf("vol_sexDAP limma.pdf", height = 6, width = 5)
  vol_sexDAP
dev.off()

seqnames <- read.csv("input files/seqnames.csv", row.names = 1)

# Retrieving the differentially bound sites
s_sex.DB <- merge(top.table_FvsM[top.table_FvsM$adj.P.Val < 0.05, ], seqnames, by = "row.names")

table(s_sex.DB$x == "groupXIX")
# FALSE  TRUE 
# 174  2091 

# in most cases, males have less accessibility than females
table(s_sex.DB$logFC > 0) # female higher accessibility if LFC is positive
# FALSE  TRUE 
# 216  2049 

# Peak at fold-difference of -1 which indicates most male sites are half as accessible as females.
hist(s_sex.DB$logFC, breaks = seq(-5, 2, 0.5))


#################### sex DAP vs DEG

# homer_annotatePeaks.sh

Peaks <- read.table("input files/ATACseq_ShiftExtend_2025_01_21_annotated_output.txt", sep = "\t", header = TRUE)
names(Peaks)[1] <- "PeakID"
DAP_sex <- Peaks[Peaks$PeakID %in% rownames(top.table_FvsM[top.table_FvsM$adj.P.Val <= 0.05, ]), ]


GE_FvsM <- read.csv("input files/DEGs_FvsM.csv", row.names = 1)
DEG_sex <- GE_FvsM[GE_FvsM$adj.P.Val < 0.05, ]

listInput <- list(
  sexDAP = unique(DAP_sex$Entrez.ID), 
  sexDEG = rownames(DEG_sex)
)

table(listInput[[1]] %in% listInput[[2]])
table(listInput[[2]] %in% listInput[[1]])

V <- Venn(listInput)
plot(V, doWeights = TRUE, doEuler = TRUE, type = "circles", show = list(Faces = FALSE))


GE_CvsE <- read.csv("input files/DEGs_CvsE.csv", row.names = 1)
DEG_CvsE <- GE_CvsE[GE_CvsE$adj.P.Val < 0.05, ]

listInput <- list(
  sexDAP = unique(DAP_sex$Entrez.ID),
  CvsEDEG = rownames(DEG_CvsE)
)

table(listInput[[1]] %in% listInput[[2]])
table(listInput[[2]] %in% listInput[[1]])


GE_Day0vs10 <- read.csv("input files/DEGs_Day0vs10.csv", row.names = 1)
DEG_trt <- GE_Day0vs10[GE_Day0vs10$adj.P.Val < 0.05, ]

listInput <- list(
  sexDAP = unique(DAP_sex$Entrez.ID),
  trtDEG = rownames(DEG_trt)
)

table(listInput[[1]] %in% listInput[[2]])
table(listInput[[2]] %in% listInput[[1]])

V <- Venn(listInput)
plot(V, doWeights = TRUE, doEuler = TRUE, type = "circles", show = list(Faces = FALSE))



#### GO

set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e111_eg58_p18/")

# sex DAP
gost.res <- gost(listInput$sexDAP[!is.na(listInput$sexDAP)], organism = "gaculeatus", user_threshold = 0.05,
                 sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
if (!is.null(gost.res)) {
  print(gostplot(gost.res))
  write.csv(gost.res$result[-c(ncol(gost.res$result))], "GO_sexDAP.csv")
} else {
  write.csv(NULL, "GO_sexDAP.csv")
}

# overlap between sex DAP and DEG
gost.res <- gost(listInput$sexDAP[listInput[[2]] %in% listInput[[1]]], organism = "gaculeatus", user_threshold = 0.05,
                 sources = c("GO:MF", "GO:CC", "GO:BP")) # correction = "fdr"
if (!is.null(gost.res)) {
  print(gostplot(gost.res))
  write.csv(gost.res$result[-c(ncol(gost.res$result))], "GO_sexDAPandDEG.csv")
} else {
  write.csv(NULL, "GO_sexDAPandDEG.csv")
}




