install.packages("knitr")
BiocManager::install("limma")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("missMethyl" ,force= TRUE, dependencies =TRUE)
BiocManager::install("DMRcate", force=TRUE)
BiocManager::install("DMRcatedata", force=TRUE, dependencies = TRUE)
#include a framework for reading in the raw data from IDAT files and 
#various specialised objects for storing and manipulating the data 
#throughout the course of an analysis
BiocManager::install("minfiData")
BiocManager::install("methylationArrayAnalysis")

#Load library
# load packages required for analysis
#we will begin with pair-wise probe analysis where each probe is tested for
#differential methylation for the comparison of interest 
library(knitr)
#limma works on a matrix data values or data that can be easily converted into 
#matrix
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(DMRcatedata)
library(methylationArrayAnalysis)
#minfi , IlluminaHumanMethylation450kanno, IlluminaHumanMethylation450kmanifest,
#missmethyl, DMRcate are methylation specific packages
destination = "D:/Bioinformatics/Nile University/Integrative/System-Biology/Methyl array/methylarray.pdf"
pdf(file=destination)
#We use limma for testing differential methylation
#2.1 Obtaining the data
# set up a path to the data directory
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
list.files(dataDirectory, recursive = TRUE)
#we downloaded the data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)
# read in the sample sheet for the experiment
#read.metharray.sheet is a minfi function that is best for reading raw meth data
#samplesheet is a csv file containing one line per sample with columns describing 
#each sample
targets <- read.metharray.sheet(dataDirectory, pattern="SampleSheet.csv")
#Basename contains the path for each IDAT file
targets
# read in the raw data from the IDAT files
#read.metharray.exp(minfi function) --> creates RGchannel Set object that belongs to
#such class
rgSet <- read.metharray.exp(targets=targets)

# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID

###########################################################################
##1st: Upstream analysis:


#2.3 Quality Control
#Now we Loaded the data and imported it into R, We can evaluate its quality
#I evaluate each CpG position in each sample based on the quality of the signal
#Using minfi package again , we can do statistical testing comparing total signal (M+U)for 
#each probe to the background signal level we use cut off p-value < 0.01
 
# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)
# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2,
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal,
       bg="white")
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group,
         pdf="qcReport.pdf")
# remove poor quality samples
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet
# remove poor quality samples from targets data
targets <- targets[keep,]
targets[,1:5]

# remove poor quality samples from detection p-value table
detP <- detP[,keep]
dim(detP)

#######################################################

#2.4 Normalization
#to minimize unwanted variations within and between samples we perform normalization
#There are different techniqyes; however, there is no universal method used ,
#In methyl arrays we use minfi framework (preprocessfunnorm) for global methylation
#such as normal/cancer
#Where preprocessQuantile used for data sets where I do not expect global diff
#betweeen samples

# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet)
# create a MethylSet object from the raw data for plotting
mSetRaw <- preprocessRaw(rgSet)
# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
#After normalization , the data is housed as GenmoicRatioObject

#######################################################################3
#2.5Exploring the data
#MDS (Multidimensional scaling) is an excellent way to view the data. THey are
#based on Principal component where each PC represent variation in the data.

# MDS plots to look at largest sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)])
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.7)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Source)])
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.7)

# Examine higher dimensions to look at other sources of variation
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(1,3))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(2,3))
legend("topleft", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], dim=c(3,4))
legend("topright", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.7, bg="white")

##########################################
#2.6 Filtering

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# exclude cross reactive probes
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
table(keep)

mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt

par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], cex=0.8)
legend("right", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.65, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Source)])
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

par(mfrow=c(1,3))
# Examine higher dimensions to look at other sources of variation
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Source)], dim=c(1,3))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Source)], dim=c(2,3))
legend("topright", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Source)], dim=c(3,4))
legend("right", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.7, bg="white")


# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values",
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values",
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))
############################################################################
#2.7 Probe-wise differential methylation analysis

# this is the factor of interest
cellType <- factor(targets$Sample_Group)
# this is the individual effect that we need to account for
individual <- factor(targets$Sample_Source)

# use the above to create a design matrix
design <- model.matrix(~0+cellType+individual, data=targets)
colnames(design) <- c(levels(cellType),levels(individual)[-1])

# fit the linear model
fit <- lmFit(mVals, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(naive-rTreg,
                            naive-act_naive,
                            rTreg-act_rTreg,
                            act_naive-act_rTreg,
                            levels=design)
contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)

write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)
# plot the top 4 most significantly differentially methylated CpGs
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=targets$Sample_Group, ylab = "Beta values")
})
##############################################################################
#2.8 Differential methylation analysis of regions
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M",
                             analysis.type = "differential", design = design,
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "naive - rTreg", arraytype = "450K")

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges
results_naive_rTreg <- as.data.frame(results.ranges[results.ranges$min_smoothed_fdr < 0.05 & results.ranges$meandiff > 0.20,])
write.csv(results_naive_rTreg,"First Condition.csv", sep=",", row.names=FALSE)
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 2, CpGs = bVals, phen.col = cols,
         what = "Beta", arraytype = "450K", genome = "hg19")

#############################################
#Performing DMR for 2nd condition
myAnnotation2 <- cpg.annotate(object = mVals, datatype = "array", what = "M",
                             analysis.type = "differential", design = design,
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "naive - act_naive", arraytype = "450K")

DMRs2 <- dmrcate(myAnnotation2, lambda=1000, C=2)
results.ranges2 <- extractRanges(DMRs)
results.ranges2
results_naive_act_naive <- as.data.frame(results.ranges2[results.ranges2$min_smoothed_fdr < 0.05 & results.ranges2$meandiff > 0.20,])
write.csv(results_naive_act_naive,"Second Condition.csv", sep=",", row.names=FALSE)
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges2, dmr = 2, CpGs = bVals, phen.col = cols,
         what = "Beta", arraytype = "450K", genome = "hg19")
######################################################################
#Condition3
myAnnotation3 <- cpg.annotate(object = mVals, datatype = "array", what = "M",
                              analysis.type = "differential", design = design,
                              contrasts = TRUE, cont.matrix = contMatrix,
                              coef = "rTreg - act_rTreg", arraytype = "450K")

DMRs3 <- dmrcate(myAnnotation3,pcutoff = 0.05, lambda=1000, C=2)
results.ranges3 <- extractRanges(DMRs3)

results_rTreg_act_rTreg <- as.data.frame(results.ranges3[results.ranges3$min_smoothed_fdr < 0.05 & results.ranges3$meandiff > 0.20,])
write.csv(results_rTreg_act_rTreg,"Third Condition.csv", sep=",", row.names=FALSE)
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges3, dmr = 2, CpGs = bVals, phen.col = cols,
         what = "Beta", arraytype = "450K", genome = "hg19")
##########################################################################################################

myAnnotation4 <- cpg.annotate(object = mVals, datatype = "array", what = "M",
                              analysis.type = "differential", design = design,
                              contrasts = TRUE, cont.matrix = contMatrix,
                              coef = "act_naive - act_rTreg", arraytype = "450K")

DMRs4 <- dmrcate(myAnnotation4, lambda=1000, C=2)
results.ranges4 <- extractRanges(DMRs4)

results_act_naive_act_rTreg <- as.data.frame(results.ranges4[results.ranges4$min_smoothed_fdr < 0.05 & results.ranges4$meandiff > 0.20,])
write.csv(results_act_naive_act_rTreg,"Fourth Condition.csv", sep=",", row.names=FALSE)
# set up the grouping variables and colours
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges4, dmr = 2, CpGs = bVals, phen.col = cols,
         what = "Beta", arraytype = "450K", genome = "hg19")





dev.off() 


