## TagSeq Limma
# last run 2025 05 27
# written by Jason Keagy

# R version 4.5.0 (2025-04-11)

# libraries
library(limma) # version 3.64.0
library(edgeR) # version 4.6.2
library(Vennerable) # version 3.1.0.9000
library(circlize) # version 0.4.16
library(ComplexHeatmap) # version 2.25.1
library(chron) # version 2.3-62
library(ggplot2) # version 3.5.2


# Read in gene counts (if one gff file for each sample)
# set working directory to folder "individual_gff_files"
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/individual_gff_files"))

file_list <- list.files(pattern = "gff")
dataset <- read.table(file_list[1], row.names = 1)
colnames(dataset)[1] <- gsub("_htseqgene_counts.gff", "", file_list[1])
for (i in 2:(length(file_list))) {     
  temp_dataset <- read.table(file_list[i], row.names = 1)
  dataset <- cbind(dataset, temp_dataset)
  colnames(dataset)[i] <- gsub("_htseqgene_counts.gff", "", file_list[i])
  rm(temp_dataset)
}
dataset <- dataset[1:(nrow(dataset) - 5), ] # only needed if used htseq

######################


# filter out the bad sequencing 
hist(colSums(dataset))
colnames(dataset[colSums(dataset) < 500])
# all six lanes of 25E6 had almost no reads
dataset <- dataset[colSums(dataset) > 500]


# Create DGEList object
d0 <- DGEList(dataset)


# Add sample information using filename
d0$samples$sampleID <- 
  gsub("S", "", gsub("-", "", gsub("_", "", substr(rownames(d0$samples), 5, 10))))
d0$samples$lane <- 
  substr(rownames(d0$samples), nchar(rownames(d0$samples)) - 7, nchar(rownames(d0$samples)) - 7)
d0$samples$ExpCon <- as.factor(substr(d0$samples$sampleID, 3, 3))
d0$samples$Fry.Number <- as.factor(substr(d0$samples$sampleID, 4, 5))


#######################
# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Add additional sample info

file1 <- "input files/Fry Brains 2018 for tagSeq.csv"
sampleInfo <- read.csv(file1)
sampleInfo$sampleID <- paste(substr(sampleInfo$Tank, 1, 2), sampleInfo$Fry.Number, sep = "")
sampleInfo$Plate <- substr(sampleInfo$Submitted.Sample.Name, 1, 1)
d0$samples <- merge(d0$samples, sampleInfo[!(colnames(sampleInfo) %in% c("Submitted.Sample.Name", "Fry.Number"))], by = "sampleID", sort = FALSE)

file2 <- "input files/Fry Brains 2018.csv"
expInfo <- read.csv(file2)
expInfo$sampleID <- paste(substr(expInfo$Tank, 1, 2), expInfo$Fry.Number, sep = "")
d0$samples <- merge(d0$samples, expInfo[c("sampleID", "Sex1", "Time.In", "Time.Attack", "Dissection.START", "Dissection.END", "Pool", "Length", "Mass")], by.x = "sampleID", sort = FALSE)

d0$samples$Sex1[d0$samples$Sex1 == "M?"] <- "M"
d0$samples$Sex1[d0$samples$Sex1 == "F?"] <- "F"
d0$samples$Sex1 <- droplevels(as.factor(d0$samples$Sex1))
d0$samples$Trt <- gsub(" ", "", as.character(d0$samples$Trt))

d0$samples$Father <- substr(d0$samples$Tank, 4, 6)

d0$samples$Dissection.START <- times(paste(as.character(d0$samples$Dissection.START), ":00", sep = ""))


# Sum technical replicates
d <- sumTechReps(d0, ID = d0$samples$sampleID)


table(d$samples$Trt, d$samples$ExpCon)
#        C  E
# Day0  64 59
# Day10 57 61


# Specify the model to be fitted. We do this before using voom since voom uses variances of the model residuals (observed - fitted)
d$samples$Group <- factor(paste0(d$samples$Trt, "_", d$samples$ExpCon, "_", d$samples$Sex1))
mm <- model.matrix(~ 0 + d$samples$Group + d$samples$SIDE.RNA) 
# mm <- model.matrix(~ d$samples$Trt * d$samples$ExpCon * d$samples$Sex1 + d$samples$SIDE.RNA)
# mm <- model.matrix(~ d$samples$Trt + d$samples$ExpCon + d$samples$Sex1 + d$samples$SIDE.RNA)
# mm <- model.matrix(~ 0 + d$samples$Group)
colnames(mm) <- gsub("Group", "", gsub("d[$]samples[$]", "", colnames(mm)))


# The next step is to remove rows that consistently have zero or very low counts. One can for example use
keep <- filterByExpr(d, mm, min.count = 10) # filter using model matrix specifications of cell group sizes
d <- d[keep, , keep.lib.sizes = FALSE]
# keep <- filterByExpr(d)
# d <- d[keep, ]
dim(d)


# Calculate normalization factors
d <- calcNormFactors(d)


# Voom
v <- voom(d, mm, plot = TRUE)
write.csv(v$E, "Limma_normalized.counts.csv")
write.csv(v$targets, "Limma_normalized.counts_sampledata.csv")


# random effect of tank:
# cor <- duplicateCorrelation(v, mm, block = d$samples$Tank)
# cor$consensus # 0.08659984
# random effect of father
# cor <- duplicateCorrelation(v, mm, block = d$samples$Father)
# cor$consensus # 0.06030443


# lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(v, mm)
# with random effect:
# fit <- lmFit(v, mm, block = d$samples$Tank, correlation = cor$consensus)
# fit <- lmFit(v, mm, block = d$samples$Father, correlation = cor$consensus)


colnames(mm)

contrasts <- makeContrasts(
  contrast_0vs10 = 
    (Day0_C_F + Day0_C_M + Day0_E_F + Day0_E_M) / 4 - 
    (Day10_C_F + Day10_C_M + Day10_E_F + Day10_E_M) / 4,
  contrast_CvsE =
    (Day0_C_F + Day0_C_M + Day10_C_F + Day10_C_M) / 4 - 
    (Day0_E_F + Day0_E_M + Day10_E_F + Day10_E_M) / 4, 
  contrast_FvsM =
    (Day0_C_F + Day0_E_F + Day10_C_F + Day10_E_F) / 4 - 
    (Day0_C_M + Day0_E_M + Day10_C_M + Day10_E_M) / 4, 
  contrast_0vs10_SexInt =
    (((Day0_C_F + Day0_E_F) / 2) - ((Day10_C_F + Day10_E_F) / 2)) -
    (((Day0_C_M + Day0_E_M) / 2) - ((Day10_C_M + Day10_E_M) / 2)),
  levels = mm)


fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)


top.table_0vs10 <- topTable(fit2, coef = "contrast_0vs10", sort.by = "P", n = Inf)
length(which(top.table_0vs10$adj.P.Val < 0.05))
write.csv(top.table_0vs10, "DEGs_Day0vs10.csv")

top.table_CvsE <- topTable(fit2, coef = "contrast_CvsE", sort.by = "P", n = Inf)
length(which(top.table_CvsE$adj.P.Val < 0.05))
write.csv(top.table_CvsE, "DEGs_CvsE.csv")

top.table_FvsM  <- topTable(fit2, coef = "contrast_FvsM", sort.by = "P", n = Inf)
length(which(top.table_FvsM$adj.P.Val < 0.05))
write.csv(top.table_FvsM, "DEGs_FvsM.csv")

top.table_0vs10_SexInt <- topTable(fit2, coef = "contrast_0vs10_SexInt", sort.by = "P", n = Inf)
length(which(top.table_0vs10_SexInt$adj.P.Val < 0.05))
write.csv(top.table_0vs10_SexInt, "DEGs_Day0vs10_SexInt.csv")

FvsM_LFC <- rownames(top.table_FvsM[top.table_FvsM$adj.P.Val < 0.05, ])
FvsM_LFC_1 <- rownames(top.table_FvsM[top.table_FvsM$logFC > 0.9 & top.table_FvsM$logFC < 1.1 &
  top.table_FvsM$adj.P.Val < 0.05, ])

stickle_gtf <- read.csv("input files/stickle.gtf", sep = "\t", header = FALSE)

#LG19 (XIX) = sex chromosome

table(stickle_gtf$V1[stickle_gtf$V9 %in% FvsM_LFC_1])

chr_FvsM <- table(stickle_gtf$V1[stickle_gtf$V9 %in% FvsM_LFC])
chr_FvsM


size = 18
vol_0vs10 <- ggplot(top.table_0vs10, aes(x = logFC, y = -log10(adj.P.Val))) +
         geom_point(data = top.table_0vs10[top.table_0vs10$adj.P.Val > 0.05, ],
            color = "grey", shape = 19, size = 2, alpha = 0.5) +
         geom_point(data = top.table_0vs10[top.table_0vs10$adj.P.Val <= 0.05 &
            top.table_0vs10$logFC > 0, ],
            color = "blue", shape = 19, size = 2, alpha = 0.5) +
         geom_point(data = top.table_0vs10[top.table_0vs10$adj.P.Val <= 0.05 &
            top.table_0vs10$logFC < 0, ],
            color = "red", shape = 19, size = 2, alpha = 0.5) +
         xlab(expression(log[2]~fold~change)) + 
         ylab(expression(-log[10]~(FDR))) +
         xlim(c(-4.5, 1.5)) + 
         ylim(c(0, 130)) + 
         geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
         theme_bw() +
         theme(axis.title = element_text(size = size),
           axis.text = element_text(size = size - 2),
           legend.title = element_text(size = size), 
           legend.text = element_text(size = size - 2))
vol_0vs10

pdf("vol_0vs10.pdf", height = 6, width = 5)
  vol_0vs10
dev.off()

size = 18
vol_FvsM <- ggplot(top.table_FvsM, aes(x = logFC, y = -log10(adj.P.Val))) +
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
  xlim(c(-4.5, 1.5)) + 
  ylim(c(0, 130)) +
         geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
         theme_bw() +
         theme(axis.title = element_text(size = size),
           axis.text = element_text(size = size - 2),
           legend.title = element_text(size = size), 
           legend.text = element_text(size = size - 2))
vol_FvsM

pdf("vol_FvsM.pdf", height = 6, width = 5)
  vol_FvsM
dev.off()


size = 18
vol_CvsE <- ggplot(top.table_CvsE, aes(x = logFC, y = -log10(adj.P.Val))) +
         geom_point(data = top.table_CvsE[top.table_CvsE$adj.P.Val > 0.05, ],
            color = "grey", shape = 19, size = 2, alpha = 0.5) +
         geom_point(data = top.table_CvsE[top.table_CvsE$adj.P.Val <= 0.05 &
            top.table_CvsE$logFC > 0, ],
            color = "blue", shape = 19, size = 2, alpha = 0.5) +
         geom_point(data = top.table_CvsE[top.table_CvsE$adj.P.Val <= 0.05 &
            top.table_CvsE$logFC < 0, ],
            color = "red", shape = 19, size = 2, alpha = 0.5) +
         xlab(expression(log[2]~fold~change)) + 
         ylab(expression(-log[10]~(FDR))) +
  xlim(c(-4.5, 1.5)) + 
  ylim(c(0, 130)) +
         geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
         theme_bw() +
         theme(axis.title = element_text(size = size),
           axis.text = element_text(size = size - 2),
           legend.title = element_text(size = size), 
           legend.text = element_text(size = size - 2))
vol_CvsE

pdf("vol_CvsE.pdf", height = 6, width = 5)
  vol_CvsE
dev.off()


listInput <- list(
  Day0vsDay10 = rownames(top.table_0vs10[top.table_0vs10$adj.P.Val < 0.05, ]),
  PoolvsBaseline = rownames(top.table_CvsE[top.table_CvsE$adj.P.Val < 0.05, ]),
  FvsM = rownames(top.table_FvsM[top.table_FvsM$adj.P.Val < 0.05, ])
)

V <- Venn(listInput)
plot(V, doWeights = TRUE, doEuler = TRUE, type = "circles", show = list(Faces = FALSE))

Sex_and_PC_intersection <- rownames(top.table_FvsM)[top.table_FvsM$adj.P.Val < 0.05] %in% rownames(top.table_0vs10)[top.table_0vs10$adj.P.Val < 0.05]

Sex_and_PC <- top.table_FvsM[top.table_FvsM$adj.P.Val < 0.05, ][Sex_and_PC_intersection, ]

Sex_and_PC <- merge(Sex_and_PC, top.table_0vs10, by = "row.names", suffixes = c(".Sex",".PC"))

plot(Sex_and_PC$logFC.Sex, Sex_and_PC$logFC.PC,
  xlab = "LFC Sex", ylab = "LFC Parental Care")
abline(a = 0, b = 1, col  = "blue")
abline(a = 0, b = -1, lty = "dashed", col = "red")
abline(v = 1, lty = "dotted")

cor.test(Sex_and_PC$logFC.Sex, Sex_and_PC$logFC.PC)
cor.test(sign(Sex_and_PC$logFC.Sex), sign(Sex_and_PC$logFC.PC))
table(sign(Sex_and_PC$logFC.Sex), sign(Sex_and_PC$logFC.PC))

chisq.test(table(sign(Sex_and_PC$logFC.Sex), sign(Sex_and_PC$logFC.PC)))


# M vs F
Expr <- scale(v$E[rownames(v$E) %in% rownames(top.table_FvsM)[which(top.table_FvsM$adj.P.Val < 0.05 & abs(top.table_FvsM$logFC) > 1)], ], scale = FALSE)

col_fun = colorRamp2(c(min(Expr), 0, max(Expr)), c("blue", "white", "red"))

ha = HeatmapAnnotation(Sex = v$targets$Sex1, CareTrt = v$targets$Trt, ExpCon = v$targets$ExpCon,
                       col = list( Sex = c("F" = "pink", "M" = "lightblue"),
                                   CareTrt = c("Day0" = "black", "Day10" = "grey90"),
                                   ExpCon = c("E" = "red", "C" = "green")))

Heatmap(Expr, col = col_fun,
        show_column_dend = FALSE,
        show_row_dend = TRUE,
        name = "Expression",
        column_title = "Samples",
        row_title = "Genes",
        show_row_names = FALSE,
        show_column_names = FALSE,
        column_split = v$targets$Group,
        top_annotation = ha,
        #width = unit(8, "cm"), height = unit(8, "cm")
)
