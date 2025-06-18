## Fry Behavior
# last run 2025 05 27
# written by Jason Keagy

# R version 4.5.0 (2025-04-11)

# libraries
library(lme4) # version 1.1-37
library(lmerTest) # version 3.1-3
library(plyr) # version 1.8.9
library(ggplot2) # version 3.5.2
library(chron) # version 2.3-62
library(car) # version 3.1-3
library(rmcorr) # version 0.7.0


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


# set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


########################### Open Field Assay

filebefore <- "input files/Fry Behavior 2018_before.csv"
resultsbefore <- read.csv(filebefore)

fileafter <- "input files/Fry Behavior 2018_after.csv"
resultsafter <- read.csv(fileafter)

results <- rbind(resultsbefore, resultsafter)

results$FishTankID <- paste(results$Tank, results$Fish.ID, sep = "_")

file3 <- "input files/Nesting Tanks 2017.csv"
fryInfo <- read.csv(file3)
fryInfo$Family <- paste(sprintf("%02d", fryInfo$Tank), sprintf("%03d", fryInfo$Male.Number), sep = "-")

file4 <- "input files/Fry Brains 2018.csv"
frySize <- read.csv(file4)
frySize$FishTankID <- paste(frySize$Tank, substr(frySize$Fry.Number, 1, 3), 
	sep = "_")
frySize$Father <- substr(frySize$Tank, 4, 6)
frySize$Sex1[frySize$Sex1 == "M?"] <- "M"
frySize$Sex1[frySize$Sex1 == "F?"] <- "F"
frySize$Sex1 <- droplevels(as.factor(frySize$Sex1))

xx <- merge(fryInfo[, c("Family", "Trt", "Clutch.Size", "Egg.Date", "Pair.Tank")], frySize[, c("Tank", "FishTankID", "Length", "Mass", "Pool", "Time.In", "Father", "SIDE.RNA", "Sex1")], by.x = "Family", by.y = "Tank")
# actual sample sizes
table(xx[!(is.na(xx$SIDE.RNA)), ]$Trt, is.na(xx[!(is.na(xx$SIDE.RNA)), ]$Pool))   
#           exp  con
# Day 0     65   65
# Day 10    62   62
table(xx[!(is.na(xx$SIDE.RNA)), ]$Trt, is.na(xx[!(is.na(xx$SIDE.RNA)), ]$Pool), 
  xx[!(is.na(xx$SIDE.RNA)), ]$Sex1)  
# , ,  = F
# 
#         
#             exp  con
#   Day 0     36   26
#   Day 10    36   29
# 
# , ,  = M
# 
#         
#             exp  con
#   Day 0     29   39
#   Day 10    26   33

frySize <- frySize[substr(frySize$Fry.Number, 1, 1) == "E", ]
mergedData <- merge(results, fryInfo[, c("Family", "Trt", "Clutch.Size", "Egg.Date", "Pair.Tank")], by.x = "Tank", by.y = "Family")
mergedData <- merge(mergedData, frySize[, c("FishTankID", "Length", "Mass", "Pool", "Time.In", "Father", "SIDE.RNA", "Sex1")], by = "FishTankID")

mergedData$Egg.Date <- as.Date(mergedData$Egg.Date, "%m/%d/%y")
mergedData$Date <- as.Date(mergedData$Date, "%m/%d/%y")
mergedData$Julian_Egg.Date <- as.numeric(format(mergedData$Egg.Date, "%j"))
mergedData$Age <- as.numeric(mergedData$Date - mergedData$Egg.Date)

# something wrong with 05-085_E2 after
mergedData <- mergedData[mergedData$FishTankID != "05-085_E2", ]

mergedData$Time.Period <- factor(mergedData$Time.Period, levels = c("before", "after"))

mergedData$Time.In <- times(paste(as.character(mergedData$Time.In), ":00", sep = ""))

mergedData$propFreeze <- mergedData$Freeze / 
	(mergedData$Freeze + mergedData$Unfreeze)
mergedData$propPlant <- mergedData$Plant / 
	(mergedData$Plant + mergedData$Open)
mergedData$propEdge <- mergedData$Edge / 
  (mergedData$Edge + mergedData$Away.From.Edge)

mergedData$Trt <- as.factor(mergedData$Trt)

mergedData$TimeBin <- factor(paste0(mergedData$Time.Period, "_", mergedData$Bin),
                        levels = c("before_0", "before_1", "after_0", "after_1"))

mergedData$Observer[mergedData$Observer == "AF "] <- "AF"
mergedData$Observer[mergedData$Observer == "BF "] <- "BF"
mergedData$Observer[mergedData$Observer == "SZ"] <- "GZ"
mergedData$Observer <- as.factor(mergedData$Observer)

mergedData$Recorder[mergedData$Recorder == " AF"] <- "AF"
mergedData$Recorder[mergedData$Recorder == " AF "] <- "AF"
mergedData$Recorder[mergedData$Recorder == "AF "] <- "AF"
mergedData$Recorder[mergedData$Recorder == " BF"] <- "BF"
mergedData$Recorder[mergedData$Recorder == " BF "] <- "BF"
mergedData$Recorder[mergedData$Recorder == "BF "] <- "BF"
mergedData$Recorder[mergedData$Recorder == ""] <- "GZ"
mergedData$Recorder[mergedData$Recorder == "JK"] <- "JCK"
mergedData$Recorder <- as.factor(mergedData$Recorder)

# Sanity check - did the experiment affect behavior?
t.test(Exploration ~ Time.Period, data = mergedData)
t.test(Exploration ~ Time.Period, data = mergedData)
t.test(propPlant ~ Time.Period, data = mergedData)
t.test(propFreeze ~ Time.Period, data = mergedData)
t.test(propEdge ~ Time.Period, data = mergedData)
t.test(Looking ~ Time.Period, data = mergedData)

# for figures
fontsize <- 24
dodge_jitter <- 0.4
dodge_stat <- 0.3
colors <- c("#53B7E0", "#E0A473")

#################### Plant ####################

# five minutes before and five minutes after
lmer.plant <- lmer(propPlant ~ Trt *  Sex1 * Time.Period + (1 | Father/Trt/FishTankID) + (1 | Observer),
               data = mergedData[!(is.na(mergedData$SIDE.RNA)) & 
                          (mergedData$TimeBin == "before_1" | mergedData$TimeBin == "after_0"), ],
               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
summary(lmer.plant)
anova(lmer.plant, type = "II")
ranova(lmer.plant)
hist(resid(lmer.plant))
shapiro.test(resid(lmer.plant))

# treatment:sex over 2 time periods
figTrtXPlant <- ggplot(data = mergedData[!(is.na(mergedData$SIDE.RNA)) & 
                          (mergedData$TimeBin == "before_1" | mergedData$TimeBin == "after_0"), ], 
  aes(x = Time.Period, y = propPlant, color = Trt, shape = Sex1)) 
figTrtXPlant <- figTrtXPlant + #scale_size_continuous(guide = "none") + 
  geom_point(aes(group = FishTankID), size = 1, alpha = 0.1, 
             position = position_dodge(width = dodge_jitter)) +
  geom_line(aes(group = FishTankID), linewidth = 1, alpha = 0.1, 
            position = position_dodge(width = dodge_jitter)) +
  stat_summary(aes(y = propPlant, group = Trt:Sex1), fun.y = mean, geom = "line",
    position = position_dodge(width = dodge_stat)) +
  stat_summary(aes(y = propPlant, group = Trt:Sex1), fun.y = mean, 
    fun.min = meanminus, fun.max = meanplus, geom = "pointrange", size = 1.25, 
    position = position_dodge(width = dodge_stat), alpha = 0.9) + 
  ylab("Proportion Time Under Plant") + 
  scale_x_discrete(labels = c("before attack", "after attack")) +
  xlab("") + 
  theme_bw() +
  ylim(c(0, 1)) +
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) +
  scale_shape_discrete(name = "Sex") +
  theme(axis.text = element_text(size = (fontsize * 0.75)),
    axis.title = element_text(size = fontsize), 
    legend.text = element_text(size = (fontsize * 0.75)),
    legend.title = element_text(size = fontsize))
figTrtXPlant

pdf("TrtXPlant.pdf", height = 6, width = 8)
  figTrtXPlant
dev.off()


# treatment:sex over all 4 time periods
figTrtXPlant <- ggplot(data = mergedData[!(is.na(mergedData$SIDE.RNA)), ], 
                        aes(x = TimeBin, y = propPlant, 
                            color = Trt, shape =Sex1)) 
figTrtXPlant <- figTrtXPlant + scale_size_continuous(guide = "none") + 
  geom_point(aes(group = FishTankID), size = 1, alpha = 0.25, 
             position = position_dodge(width = dodge_jitter)) +
  geom_line(aes(group = FishTankID), linewidth = 1, alpha = 0.1, 
            position = position_dodge(width = dodge_jitter)) +
  stat_summary(aes(group = Trt:Sex1), fun.y = mean, geom = "line",
               position = position_dodge(width = dodge_stat)) +
  stat_summary(aes(y = propPlant, group = Trt:Sex1), fun.y = mean, 
    fun.min = meanminus, fun.max = meanplus, geom = "pointrange", size = 1.25, 
    position = position_dodge(width = dodge_stat), alpha = 0.9) +
  ylab("Proportion Time Under Plant in\nOpen Field Assay") +
  xlab("Time Period (min since entered pool)") + 
  scale_x_discrete(labels = c("0-5", "5-10", "10-15", "15-20")) + 
  theme_bw() + ylim(c(0, 1)) + 
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) +
  scale_shape_discrete(name = "Sex") +
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize))
figTrtXPlant


########################### Social Experiment

file1 <- "input files/Social Experiment position data.csv"
dist <- read.csv(file1, stringsAsFactors = TRUE)
dist$FryID <- paste(dist$Tank, dist$Fry.Number, sep = "-")

file2 <- "input files/Nesting Tanks 2017.csv"
fryInfo <- read.csv(file2, stringsAsFactors = TRUE)
fryInfo$Family <- fryInfo$Family <- paste(sprintf("%02d", fryInfo$Tank), sprintf("%03d", fryInfo$Male.Number), sep = "-")

file3 <- "input files/Social Experiment.csv"
expInfo <- read.csv(file3, stringsAsFactors = TRUE)
expInfo$FryID <- paste(expInfo$Tank, expInfo$Fry.Number, sep = "-")

mergedData <- merge(dist, fryInfo[, c("Family", "Male.Number", "Trt", "Clutch.Size", "Egg.Date", "Pair.Tank")], by.x = "Tank", by.y = "Family")
mergedData <- merge(mergedData, expInfo[, c("FryID", "Mass", "Length", "Density", "Sex")], by = "FryID")

mergedData$Time.Period <- relevel(mergedData$Time.Period, ref = "Before")

# get observer ID
mergedData$Obs <- gsub("\\.", "", substr(mergedData$File, 10, 12))

# one of the photos probably had the wrong scale because it has impossible distances
mergedData <- mergedData[mergedData$Dist.Edge < 700, ]
mergedData <- mergedData[mergedData$Dist.Moved < 10000 | is.na(mergedData$Dist.Moved), ]

# new dataset with means for each fish for each period
DistFish <- ddply(mergedData, c("FryID", "Male.Number", "Tank", "Time.Period", "Trt", "Sex", "Obs"), 
                  summarise, meanMinDistFish = mean(Min.Dist.Fish, na.rm = TRUE))

hist(DistFish$meanMinDistFish)
hist(log(DistFish$meanMinDistFish))
hist(sqrt(DistFish$meanMinDistFish))
shapiro.test(DistFish$meanMinDistFish)
shapiro.test(log(DistFish$meanMinDistFish))
shapiro.test(sqrt(DistFish$meanMinDistFish))

DistFish$ln.meanMinDistFish <- log(DistFish$meanMinDistFish)
leveneTest(ln.meanMinDistFish ~ Trt * Sex * Time.Period, data = DistFish)

l <- lmer(ln.meanMinDistFish ~ Trt * Sex * Time.Period + (1 | Male.Number/Tank/FryID) + (1 | Obs), 
          data = DistFish)
summary(l)
anova(l, type = "II")
ranova(l)
hist(residuals(l))
shapiro.test(residuals(l))


fig_MinDistFish <- ggplot(DistFish, aes(Time.Period, meanMinDistFish, 
                                        color = Trt, shape = Sex))
fig_MinDistFish <- fig_MinDistFish + 
  geom_point(aes(group = FryID), size = 1, alpha = 0.25, 
             position = position_dodge(width = dodge_jitter)) +
  geom_line(aes(group = FryID), linewidth = 1, alpha = 0.1, 
            position = position_dodge(width = dodge_jitter)) +
  stat_summary(aes(group = Trt:Sex), fun = mean, geom = "line",
               position = position_dodge(width = dodge_stat)) +
  stat_summary(aes(y = meanMinDistFish, group = Trt:Sex), fun.y = mean, 
               fun.ymin = meanminus, fun.ymax = meanplus, geom = "pointrange", size = 1.25, 
               position = position_dodge(width = dodge_stat), alpha = 0.9) +
  theme_bw() +
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) +
  scale_shape_discrete(name = "Sex") +
  scale_x_discrete(labels = c("before attack", "after attack")) +
  scale_y_continuous(transform = "log", breaks = c(0, 25, 50, 75, 100, 125)) +
  xlab("") + 
  ylab("Nearest Fish Distance (mm)") +
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize))
fig_MinDistFish

pdf("MinDistFish.pdf", height = 6, width = 8)
  fig_MinDistFish
dev.off()


# new dataset with means for each fish for each period
DistMoved <- ddply(mergedData, c("FryID", "Male.Number", "Tank", "Time.Period", "Trt", "Sex", "Obs"), 
                   summarise, meanDistMoved = mean(Dist.Moved, na.rm = T))

hist(DistMoved$meanDistMoved)
hist(log(DistMoved$meanDistMoved))
hist(sqrt(DistMoved$meanDistMoved))
shapiro.test(DistMoved$meanDistMoved)
shapiro.test(log(DistMoved$meanDistMoved))
shapiro.test(sqrt(DistMoved$meanDistMoved))

DistMoved$sqrt.meanDistMoved <- sqrt(DistMoved$meanDistMoved)
leveneTest(sqrt.meanDistMoved ~ Trt * Sex * Time.Period, data = DistMoved)

l <- lmer(sqrt.meanDistMoved ~ Trt * Sex * Time.Period + (1 | Male.Number/Tank/FryID) + (1 | Obs), 
          data = DistMoved, 
          control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
summary(l)
anova(l, type = "II")
ranova(l)


fig_DistMoved <- ggplot(DistMoved, aes(Time.Period, meanDistMoved, 
                                       color = Trt, shape = Sex))
fig_DistMoved <- fig_DistMoved + 
  geom_point(aes(group = FryID), size = 1, alpha = 0.25, 
             position = position_dodge(width = dodge_jitter)) +
  geom_line(aes(group = FryID), linewidth = 1, alpha = 0.1, 
            position = position_dodge(width = dodge_jitter)) +
  stat_summary(aes(group = Trt:Sex), fun.y = mean, geom = "line",
               position = position_dodge(width = dodge_stat)) +
  stat_summary(aes(y = meanDistMoved, group = Trt:Sex), fun.y = mean, 
               fun.ymin = meanminus, fun.ymax = meanplus, geom = "pointrange", size = 1.25, 
               position = position_dodge(width = dodge_stat), alpha = 0.9) +
  theme_bw() +
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) +
  scale_shape_discrete(name = "Sex") +
  scale_x_discrete(labels = c("before attack", "after attack")) +
  scale_y_continuous(transform = "sqrt", breaks = c(0, 50, 100, 150, 200, 250)) +
  xlab("") + 
  ylab("Activity (mm/10s)") +
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize))
fig_DistMoved

pdf("DistMoved.pdf", height = 6, width = 8)
  fig_DistMoved
dev.off()


########################### Scototaxis

file <- "input files/Scototaxis.csv"
scot <- read.csv(file, stringsAsFactors = TRUE)

# data manipulation

# just use paired tanks
scot$Father <- substr(scot$Tank, 4, 6)
pair_father <- names(table(scot$Father)[table(scot$Father) >= 2 * min(table(scot$Tank))])

scot_pair <- scot[scot$Father %in% pair_father, ]
scot_pair$Tank <- droplevels(scot_pair$Tank)

# if sex not determined by gonads, use sex from morphology
scot_pair$Gonad.Sex[scot_pair$Gonad.Sex == "?"] <- scot_pair$Sex[scot_pair$Gonad.Sex == "?"]
scot_pair$Gonad.Sex[is.na(scot_pair$Gonad.Sex)] <- "F"
scot_pair$Gonad.Sex <- droplevels(scot_pair$Gonad.Sex)

# format "Time in White"
scot_pair$Time.in.White <- as.numeric(as.difftime(as.character(scot_pair$Time.in.White), format = "%M:%S", units = "secs"))

# format "Time to Enter White"
scot_pair$Time.to.Enter.White <- as.numeric(as.difftime(as.character(scot_pair$Time.to.Enter.White), format = "%M:%S", units = "secs"))

# sample sizes
table(scot_pair$Tank)
table(scot_pair$Treatment:scot_pair$Gonad.Sex)
table(scot_pair$Treatment)


TinW <- lmer(Time.in.White ~ Treatment * Gonad.Sex + (1 | Father/Treatment) + (1 | Observer), 
             data = scot_pair)
# singular fit caused by random effects having no variance
summary(TinW)
anova(TinW, type = "II")

NS <- lmer(Number.Switches ~ Treatment * Gonad.Sex + (1 | Father/Treatment) + (1 | Observer), 
           data = scot_pair)
# singular fit caused by random effects having no variance
summary(NS)
anova(NS, type = "II")

TEW <- lmer(log(1 + Time.to.Enter.White) ~ Treatment * Gonad.Sex + (1 | Father/Treatment) + (1 | Observer), 
            data = scot_pair)
# singular fit caused by random effects having no variance
summary(TEW)
anova(TEW, type = "II")

# correlation table

# PCA
scot_pair$lnTime.to.Enter.White <- log(1 + scot_pair$Time.to.Enter.White)
trait_matrix <- scot_pair[c("Time.in.White", "Number.Switches", "lnTime.to.Enter.White")]
ScotPCA <- prcomp(trait_matrix, scale. = TRUE, center = TRUE)
ScotPCA
summary(ScotPCA)
scot_pair$PC1 <- ScotPCA$x[, 1]
scot_pair$PC2 <- ScotPCA$x[, 2]

hist(scot_pair$PC1)
summary(scot_pair$PC1)
shapiro.test(scot_pair$PC1)

hist(scot_pair$PC2)
summary(scot_pair$PC2)
shapiro.test(scot_pair$PC2)

PCA <- lmer(PC1 ~ Treatment * Gonad.Sex + 
              (1 | Father/Treatment) + (1 | Observer), data = scot_pair)
summary(PCA)
anova(PCA, type = "II")
# singular fit caused by random effects having no variance
hist(residuals(PCA))
shapiro.test(residuals(PCA))


# plot
PC1Plot <- ggplot(scot_pair, aes(Treatment:Gonad.Sex, PC1, color = Treatment, shape = Gonad.Sex))
PC1Plot <- PC1Plot + 
  geom_jitter(width = dodge_jitter, alpha = 0.25, size = 1) +
  stat_summary(aes(y = PC1, group = Treatment:Gonad.Sex), fun.y = mean, 
               fun.ymin = meanminus, fun.ymax = meanplus, geom = "pointrange", size = 1.25, 
               position = position_dodge(width = dodge_stat), alpha = 0.9) +
  scale_y_continuous(name = "PC1 (higher = more cautious)") +
  theme_bw() + 
  scale_color_discrete(type = colors, name = "Care\nTreatment", labels = c("Orphaned", "Father-reared")) +
  scale_x_discrete(breaks = NULL) +
  scale_shape_discrete(name = "Sex") +
  xlab("") + 
  theme(axis.text = element_text(size = (fontsize * 0.75)),
        axis.title = element_text(size = fontsize), 
        legend.text = element_text(size = (fontsize * 0.75)),
        legend.title = element_text(size = fontsize))
PC1Plot

pdf("Scototaxis PC1.pdf", height = 6, width = 5)
  PC1Plot
dev.off()

