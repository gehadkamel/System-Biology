#load packages
BiocManager::install("DESeq2", force = TRUE)#one time
library(DESeq2)#each time you open your R
library(ggplot2)
library(RColorBrewer)
library(NMF)
####################
#1 Data cleaning and exploration
#loading the gene expression data 
data = as.matrix(read.csv("DEseq2/bladder_counts.csv", row.names=1, header = T))
#loading the phenotype data 
pheno = read.csv("DEseq2/bladder_sample_sheet.csv", row.names=1)

table(pheno$Sample.Type)
#expolore data distribution using histogram
hist(data, col = "orange", main = "Histogram")
#get the log2 of the data for better visualization
hist(log2(data +1), col = "orange", main= "Histogram")
#Visualize data using boxplot
boxplot(log2(data[1:5,]+1))

#QQ plot for normality check
#if data are not aligned to my line so they are not normally distributed
qqnorm(data[30,])
qqline(data[30,])
#Explore if there is any missing expression value
sum(is.na(data))
is.na(data)
is.null(data)
is.nan(data)
#We also can use shapiro testing for normality check
####
#make sure that the columns of my data match the rows of the pheno data
pehno3= pheno[colnames(data),]

#Because DEseq data requires the data type to be only integers
#We will keep the gene names in a variable
genes = rownames(data)
#Now we can convert data into integers
data = apply(round(data), 2 , as.integer)
#Rename the rows of data to gene names
rownames(data) = genes
################################################################################
#2 Differential expression

#specify how many conditions do we have based on the column of the phenotypic table
cond1 = "Primary Tumor"
cond2= "Solid Tissue Normal"

#Create DEseq data set object
dds= DESeqDataSetFromMatrix(countData = data, colData = pheno, design = ~Sample.Type)
#Run DESeq workflow
dds.run = DESeq(dds)
#Specifying the contrast so that based on the twoconditions I have
res=results(dds.run, contrast = c("Sample.Type", cond1, cond2))
#remove nulls
# remove nulls
res=as.data.frame(res[complete.cases(res), ])


#Choose the statistical significant DEGs based on the adjusted p-value <0.05
#and log2 fold change more than 2 
deseq.deg= res[res$padj < 0.05 & abs(res$log2FoldChange)>1.5,]

#Export degs into my folder for further analysis
write.csv(as.matrix(deseq.deg), file= "DESeq DEGs.csv", quote=F, row.names = T)
#We can visualize the most  significant genes through gprofiler GUI for 
#enrichment analysis

#Draw degs volcano plot
path.VCplot <- ("D:/Bioinformatics/Nile University/Integrative/System-Biology/VC.png")
png(path.VCplot, width = 1000, height = 1000, res = 200)
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Control vs infection DEGs"))
with(subset(res, padj<.05 & (log2FoldChange)>1.5), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & (log2FoldChange)< -1.5), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-7,y=45,c("upregulated","downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=19)

dev.off()

#####################################################

###Draw heatmap#####
path.HMplot <- ("D:/Bioinformatics/Nile University/Integrative/System-Biology/HM.png")
png(path.HMplot, width= 1000, height = 1000, res =200)
#normalize the data
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))
#extract counts values of DEGs only for each stage
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
aheatmap(log2(exp.degs+1), annCol =pheno$Sample.Type, col = rev(brewer.pal(9,"RdBu")), main="mRNA Control vs infection")

heatmap(log2(exp.degs+1))


dev.off()
