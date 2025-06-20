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


# QC stuff
pdf(("Fragment_Size_Distribution_2024.pdf"))
for (i in list.files("bam_files/", pattern = "bam$", full.names = TRUE)) {
	bamfile <- i
    bamfile.labels <- gsub("_MTdups_removed_mapped.bam", "", basename(i))
    ## generate fragement size distribution
    fragSize <- fragSizeDist(bamfile, bamfile.labels)
}
dev.off()


# now the diffBind stuff
s <- dba(sampleSheet = "input files/samples_ShiftExtend.csv")

# count
ss <- dba.count(s, bUseSummarizeOverlaps = TRUE) # this can take a while!
ss

# ss2 <- dba.count(s, summits = 75) # this can take a while!
# ss2

x <- dba.show(ss)$FRiP
mean(x)
sd(x)

counts <- dba.peakset(ss, bRetrieve = TRUE)


ss_Raw <- dba.count(ss, peaks = NULL, score = DBA_SCORE_READS)
counts_Raw <- dba.peakset(ss_Raw, bRetrieve = TRUE)
dba.peakset(ss_Raw, bRetrieve = TRUE, writeFile = "AccessibilityMatrix.csv")


plot(ss, density.info = "none", cexRow = 1, cexCol = 1)
dba.plotPCA(ss, attributes = c(DBA_CONDITION, DBA_FACTOR))


# Differential binding analyis
s_patcare <- dba.contrast(ss, categories = DBA_CONDITION)
s_sex <- dba.contrast(ss, categories = DBA_FACTOR)

s_sex <- dba.analyze(s_sex, method = DBA_ALL_METHODS)
dba.show(s_sex, bContrasts = T)

s_patcare <- dba.analyze(s_patcare, method = DBA_ALL_METHODS)
dba.show(s_patcare, bContrasts = T)

plot(s_sex, contrast = 1, density.info = "none", cexRow = 1, cexCol = 1)
dba.plotHeatmap(s_sex, contrast = 1, correlations = FALSE, scale = "row")

# Retrieving the differentially bound sites
s_sex.DB <- dba.report(s_sex)

table(seqnames(s_sex.DB) == "groupXIX")
# FALSE  TRUE 
# 232  2069

# in most cases, males have less accessibility than females
table((s_sex.DB$Conc_M - s_sex.DB$Conc_F) > 0)
# FALSE  TRUE 
# 1986   315

# Peak at fold-difference of -1 which indicates most male sites are half as accessible as females.
# no dosage compensation in stickleback? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4635650/
hist(s_sex.DB$Fold, breaks = seq(-5, 2, 0.5))

# Volcano plot
dba.plotVolcano(s_sex)

sexDAP <- s_sex$contrasts[[1]]$DESeq2$de

size = 18
vol_sexDAP <- ggplot(sexDAP, aes(x = fold, y = -log10(padj))) +
  geom_point(data = sexDAP[sexDAP$padj > 0.05, ],
             color = "grey", shape = 19, size = 2, alpha = 0.5) +
  geom_point(data = sexDAP[sexDAP$padj <= 0.05 &
                             sexDAP$fold > 0, ],
             color = "blue", shape = 19, size = 2, alpha = 0.5) +
  geom_point(data = sexDAP[sexDAP$padj <= 0.05 &
                             sexDAP$fold < 0, ],
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

pdf("vol_sexDAP diffbind.pdf", height = 6, width = 5)
  vol_sexDAP
dev.off()


x <- levels(seqnames(s_sex.DB))[grepl("^group", levels(seqnames(s_sex.DB)))]
s_sex_GR <- GenomicRanges::GRanges(s_sex.DB)
# I can't figure out how to override covplot resorting everything stupidly
p <- covplot(s_sex_GR, chrs = x) + ggtitle("ATACseq Sex Differential Accessibility Peaks")
print(p)


# save dba object
dba.save(ss,'ATACseq_ShiftExtend_2025_01_21')