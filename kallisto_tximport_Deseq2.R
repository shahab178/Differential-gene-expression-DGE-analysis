#Install bioconductor 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
#Install the libraries
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicAlignments")
BiocManager::install("Rhtslib")
BiocManager::install("lasso2")
BiocManager::install("DEGreport")
BiocManager::install("PCAtools")
BiocManager::install("tximport")
BiocManager::install("tximportData")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("rhdf5")
BiocManager::install("apeglm")
BiocManager::install("DEGreport")
BiocManager::install("vsn")
# The easiest way to get readr is to install the whole tidyverse:
install.packages("tidyverse")
install.packages("readr")
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("openxlsx")
install.packages("writexl")
install.packages("ashr")
install.packages("pheatmap")
devtools::install_github("KlugerLab/DAseq")

#Load packages
library("EnsDb.Hsapiens.v86")
library("tximportData")
library("tximport")
library("readr")
library("ggplot2")
library("dplyr")
library("DESeq2")
library("plyr")
library("reshape2")
library("tidyverse")
library("biomaRt")
library("RColorBrewer") 
library("rhdf5")
library("apeglm")
library("openxlsx")
library("writexl")
library("DEGreport")
library("ashr")
library("EDASeq")
library("vsn")
library(pheatmap")

################
#liver-Syl_Vac# 
################
getwd()

rm(list= ls()[!(ls() %in% c('ensembl'))])
#rm(list = ls())
dev.off(dev.list()["RStudioGD"])

setwd("/.../liver/")

#dir to the path which everything is there
dir <- "/.../liver/results"

#dir to path where the kallisto run are stored in separated folders eg. F06, F20, etc)
sample_id <- (file.path(dir, "results"))
sample_id

#make a folder with the name of sample_table and put the sample table.csv there
sample_table <- (file.path(dir, "sample_table"))
sample_table
#If ensemble problem with connection happen you should change the miror at the end
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

txdb <- EnsDb.Hsapiens.v86
tx2gene <- transcripts(txdb, 
                       columns = c(listColumns(txdb, "tx_id"), "tx_id", "gene_id", "gene_name"),
                       return.type = "DataFrame")
head(tx2gene)
######################
sample_info_Syl_Vac <- read.csv(file.path(sample_table, "sample_table_S_V.csv"), header = TRUE, sep = ",") 
head(sample_info_Syl_Vac)
#Here we import from abundance.h5 instead of tsv file
abundance_samples_counts_Syl_Vac <-(file.path(dir,sample_info_Syl_Vac$Run,"abundance.h5"))

names(abundance_samples_counts_Syl_Vac)<-sample_info_Syl_Vac$Run
abundance_samples_counts_Syl_Vac
file.exists(abundance_samples_counts_Syl_Vac)

#set txOut=F to aggregate transcript abundances to the gene level.
txi_kallisto_Syl_Vac<-tximport(abundance_samples_counts_Syl_Vac, type = "kallisto", tx2gene = tx2gene, txOut = F, ignoreTxVersion = T)
head(txi_kallisto_Syl_Vac$counts)
names(txi_kallisto_Syl_Vac)

total_counts_samples_Syl_Vac<-txi_kallisto_Syl_Vac$counts
total_counts_samples_Syl_Vac<-as.data.frame(total_counts_samples_Syl_Vac)
nrow(total_counts_samples_Syl_Vac)
total_counts_samples_Syl_Vac<-total_counts_samples_Syl_Vac[ rowSums(total_counts_samples_Syl_Vac)!=0, ] 
nrow(total_counts_samples_Syl_Vac)

total_counts_samples_Syl_Vac<- rownames_to_column(total_counts_samples_Syl_Vac)

total_counts_samples_Syl_Vac<- as.data.frame(total_counts_samples_Syl_Vac) %>% 
  tibble::rownames_to_column("ensembl_gene_id")

liver_table_Syl_Vac <- write.table( total_counts_samples_Syl_Vac, "/.../liver/liver_Syl_Vac.txt")

######## Deseq2########Differential gene expression analysis

dds_Syl_Vac <- DESeqDataSetFromTximport(txi_kallisto_Syl_Vac, sample_info_Syl_Vac, ~treatment)
dds_Syl_Vac = estimateSizeFactors(dds_Syl_Vac)
dds_Syl_Vac = estimateDispersions(dds_Syl_Vac)
plotDispEsts(dds_Syl_Vac)
dds_Syl_Vac$treatment <- relevel(dds_Syl_Vac$treatment, ref = "Vaccine")
#for normalization of LogFC here we use "apeglm" package
dds_Syl_Vac <- DESeq(dds_Syl_Vac)
res_Syl_Vac <- lfcShrink(dds_Syl_Vac, coef = 2, type = "apeglm")

plotCounts(dds_Syl_Vac, gene = which.min(res_Syl_Vac$padj), intgroup = "treatment")
#For saving Fragments Per Kilobase per Million mapped fragments (FKPM)
fkpm <- fpkm(dds_Syl_Vac, robust = TRUE)
liver_DESeq_Syl_Vac <- write.table(as.data.frame(fkpm),"/.../liver/fkpm.txt")
#Scree plot

library('PCAtools')
vst <- assay(vst(dds_Syl_Vac))
p<- pca(vst, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
#A pairs plot
pairsplot(p)
#
elbow <- findElbowPoint(p$variance)
elbow
#A bi-plot
biplot(p)
biplot(p, showLoadings = T,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
#https://bioconductor.riken.jp/packages/3.6/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
###lfshrink test###
resultsNames(dds_Syl_Vac)
resLFC <- lfcShrink(dds_Syl_Vac, coef=2)
resLFC
############
res <- results(dds_Syl_Vac)
library("BiocParallel")
#register(MulticoreParam(4))
register(SnowParam(4))
resOrdered <- res_Syl_Vac[order(res_Syl_Vac$pvalue),]
summary(res_Syl_Vac)
#How many adjusted p-values were less than 0.1?
sum(res_Syl_Vac$padj < 0.1, na.rm=TRUE)
res05 <- results(dds_Syl_Vac, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
#BiocManager::install("IHW")
library("IHW")
resIHW <- results(dds_Syl_Vac, filterFun=ihw)
summary(resIHW) 
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult
#Exploring and exporting results
plotMA(res_Syl_Vac, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))
#Alternative shrinkage estimators

resApe <- lfcShrink(dds_Syl_Vac, coef=2, type="apeglm")
resAsh <- lfcShrink(dds_Syl_Vac, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="normal")
plotMA(resApe, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
#Remove genes that have almost no information in any of the given samples
#Optional
dds <- dds_Syl_Vac[rowSums(counts(dds_Syl_Vac)) >1,]
dds <- DESeq(dds_Syl_Vac)

##P-value Distribution
#https://isglobal-brge.github.io/Master_Bioinformatics/rnaseq-data-analysis.html#differential-expression-analysis
DEresults <- results (dds_Syl_Vac, contrast = c("treatment", 'Sylvatic', 'Vaccine'))
DEresults <- DEresults[order(DEresults$pvalue),]
DEresults

DESeq2::plotMA(object = dds_Syl_Vac, ylim = c(-5,5))
DESeq2::plotMA(object = res_Syl_Vac, ylim = c(-5,5))

ggplot(data= as.data.frame(DEresults), aes(x = pvalue)) +
  geom_histogram(bins = 100)

#RLE(relative Log Expression) plots on the raw counts and normalized counts


par(mfrow = c(2, 1))
plotRLE(DESeq2::counts(dds_Syl_Vac, normalized = FALSE), 
        outline=FALSE, ylim=c(-2, 2), 
        col = as.numeric(sample_info_Syl_Vac$treatment), 
        main = 'Raw counts')
plotRLE(DESeq2::counts(dds_Syl_Vac, normalized = TRUE), 
        outline=FALSE, ylim=c(-2, 2), 
        col = as.numeric(sample_info_Syl_Vac$treatment), 
        main = 'Normalized Counts (DESeq2)')


#obtain regularized log-transformed values
DESeq.rlog <- rlog(dds_Syl_Vac, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)
#maen-sd plot for log-transformed data

msd_Plot <- meanSdPlot(rlog.norm.counts,
                       ranks = FALSE,
                       plot = FALSE)
msd_Plot$gg + 
  ggtitle("rlog-transformed read counts") + 
  ylab("standard deviation")

#Hierarchical clustering
distance.m_rlog <- as.dist(1-cor(rlog.norm.counts, method = "pearson"))

plot(hclust(distance.m_rlog),
labels = colnames(rlog.norm.counts),
main = "rlog transformed read counts\ndistance: Pearson correlation")


rld = rlog(dds_Syl_Vac)
z <- plotPCA(rld, intgroup = c('treatment'))
z + geom_label(aes(label = name))


DGE.results <- results(dds_Syl_Vac,
                       independentFiltering = TRUE,
                       alpha = 0.05)
table(DGE.results$padj < 0.05)

#Volcano from Deseq2 
EnhancedVolcano(res_Syl_Vac,
                lab = rownames(res_Syl_Vac),
                x = 'log2FoldChange',
                y = 'pvalue')
#############################################################
resOrdered <- res_Syl_Vac[order(res_Syl_Vac$padj),]
resSig <- subset(resOrdered, padj < 0.05)
resl2fc <- resSig[order(resSig$log2FoldChange),]
resl2fc_UP <- subset(resl2fc, log2FoldChange > 1.0)
resl2fc_DOWN <- subset(resl2fc, log2FoldChange < -1.0)

#QC analysis
#https://bioinformatics-core-shared-training.github.io/Bulk_RNAseq_Course_Apr22/Bulk_RNAseq_Course_Base/Markdowns/07_Data_Exploration.html
rawCounts <- round(txi_kallisto_Syl_Vac$counts, 0) 
# check dimension of count matrix
dim(rawCounts)
# for each gene, compute total count and compare to threshold
# keeping outcome in vector of 'logicals' (ie TRUE or FALSE, or NA)
keep <- rowSums(rawCounts) > 5
# summary of test outcome: number of genes in each class:
table(keep, useNA="always")

# subset genes where test was TRUE
filtCounts <- rawCounts[keep,]
# check dimension of new count matrix
dim(filtCounts)
#Raw counts
summary(filtCounts)

# few outliers affect distribution visualization
boxplot(filtCounts, main='Raw counts', las=2)
#https://sbc.shef.ac.uk/workshops/rnaseq-r-online_v1/session1.nb.html#Quality_control_of_the_imported_counts
dds_Syl_Vac$Run
head(txi_kallisto_Syl_Vac$counts)

# Raw counts mean expression Vs standard Deviation (SD)
plot(rowMeans(filtCounts), rowSds(filtCounts), 
     main='Raw counts: sd vs mean', 
     xlim=c(0,10000),
     ylim=c(0,5000))

# Get log2 counts
logcounts <- log2(filtCounts + 1)
# summary(logcounts[,1]) # summary for first column
# summary(logcounts) # summary for each column

# make a colour vector
statusCols <- str_replace_all(sample_info_Syl_Vac$treatment, c(Vaccine="red", Sylvatic="orange"))

# Check distributions of samples using boxplots
boxplot(logcounts,
        xlab="",
        ylab="Log2(Counts)",
        las=2,
        col=statusCols,
        main="Log2(Counts)")
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(logcounts), col="blue")

# Log2 counts standard deviation (sd) vs mean expression
plot(rowMeans(logcounts), rowSds(logcounts), 
     main='Log2 Counts: sd vs mean')
#VST : variance stabilizing transformation
#Variance stabilizing transformation (VST) aims at generating a matrix of values for which variance is constant across the range of mean values, especially for low mean.
#The vst function computes the fitted dispersion-mean relation, derives the transformation to apply and accounts for library size.
vst_counts <- vst(filtCounts)

# Check distributions of samples using boxplots
boxplot(vst_counts, 
        xlab="", 
        ylab="VST counts",
        las=2,
        col=statusCols)
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(vst_counts), col="blue")

# VST counts standard deviation (sd) vs mean expression
plot(rowMeans(vst_counts), rowSds(vst_counts), 
     main='VST counts: sd vs mean')

#Principal Component Analysis
library(ggfortify)

rlogcounts <- rlog(filtCounts)

# run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
autoplot(pcDat)

autoplot(pcDat,
         data = sample_info_Syl_Vac, 
         colour="treatment", 
         shape="gender_s",
         size=5)

# setting shape to FALSE causes the plot to default to using the labels instead of points
library('ggrepel')
autoplot(pcDat,
         data = sample_info_Syl_Vac, 
         colour="treatment", 
         shape="gender_s",
         size=5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=Run), box.padding = 0.8)

sampleinfo <- mutate(sample_info_Syl_Vac, Status=case_when(
  Run=="F20" ~ "Sylvatic",
  Run=="F37" ~ "Vaccine", 
  TRUE ~ treatment))
#and export it so that we have the correct version for later use
write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")

autoplot(pcDat,
         data = sampleinfo, 
         colour="treatment", 
         shape="gender_s",
         size=5)
#total number of reads for that sample
sum(assay(dds_Syl_Vac)[,1])
#for getting reads in all samples
colSums(assay(dds_Syl_Vac))
#Visualising count distributions
boxplot(assay(dds_Syl_Vac))
#Get log2 counts
vsd <- vst(dds_Syl_Vac,blind=TRUE)
# Check distributions of samples using boxplots
boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(assay(vsd)), col="blue")
#PCA
plotPCA(vsd,intgroup="treatment")

library('dplyr')
library('ggplot2')
pca_data <- plotPCA(vsd,intgroup="treatment",returnData = TRUE) %>% 
  rename(run = name) %>% 
  left_join(sample_info_Syl_Vac)


#Import some of the result's tables
liver_DESeq_Syl_Vac <- write.table(as.data.frame(res),"/.../liver/res.txt")
liver_DESeq_Syl_Vac <- write.table(as.data.frame(res_Syl_Vac),"/.../liver/liver_DESeq_Syl_Vac.txt")
liver_DESeq_Syl_Vac_degs <- write.table(as.data.frame(resSig),"/.../liver/liver_DESeq_Syl_Vac_DEGS.txt")
liver_DESeq_Syl_Vac_up <- write.table(as.data.frame (resl2fc_UP), "/.../liver/liver_DESeq_Syl_Vac_up.txt", row.names = TRUE)
liver_DESeq_Syl_Vac_down <- write.table(as.data.frame (resl2fc_DOWN), "/.../liver/liver_DESeq_Syl_Vac_donw.txt", row.names = TRUE)


#Volcano Plot (with apelgm adjustment)
topT <- as.data.frame(res_Syl_Vac)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="res_Syl_Vac", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)))

#Volcano Plot (without apelgm adjustment)
topT <- as.data.frame(res)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="res", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))

abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)))



resBigFC <- results(dds_Syl_Vac, lfcThreshold=1, altHypothesis="greaterAbs")
plotMA(resBigFC, ylim=c(-5,5))
abline(h=c(-1,1),lwd=5)

##biomart for converting the gene's symbols

#DEGS#change the header and put ID at the first header column for up and down and DEGs
getwd()

table_DEGS=read.csv("liver_DESeq_Syl_Vac_DEGS.csv", sep = ",", header = T)
transcript_ids <- table_DEGS$ID
transcript_ids
res <- getBM(attributes = c('ensembl_gene_id',
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_gene_id',
             values = transcript_ids,
             mart = ensembl) #alterar de acordo com o objeto do biomart
res
c <- merge(table_DEGS, res, by.x = "ID", by.y = "ensembl_gene_id")
d <- c[order(c$padj),]
e <- d[!duplicated(d$external_gene_name),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]

write.xlsx(as.data.frame(c),"transcripts_converted_fkpm_liver.xlsx", rowNames=TRUE)
write.xlsx(as.data.frame(g),"transcripts_converted_DEGS_liver.xlsx", rowNames=TRUE)
write.xlsx(as.data.frame(e),"e.xlsx", rowNames=TRUE)

#For volcano with everything from DESEq2
table_DEGS=read.csv("liver_DESeq_Syl_Vac.csv", sep = ",", header = T)
transcript_ids <- table_DEGS$ID
transcript_ids
res <- getBM(attributes = c('ensembl_gene_id',
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_gene_id',
             values = transcript_ids,
             mart = ensembl) #alterar de acordo com o objeto do biomart
res
c <- merge(table_DEGS, res, by.x = "ID", by.y = "ensembl_gene_id")
d <- c[order(c$padj),]
e <- d[!duplicated(d$external_gene_name),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]

write.xlsx(as.data.frame(g),"transcripts_converted_DEGS_liver.xlsx", rowNames=TRUE)

#Some more Qc analysis
ntd <- normTransform(dds_Syl_Vac)
library("vsn")
meanSdPlot(assay(ntd))
#####
library("pheatmap")
select <- order(rowMeans(counts(dds_Syl_Vac,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_Syl_Vac)[,c("treatment", "Run")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
###
vsd <- vst(dds_Syl_Vac, blind = FALSE)
rld <- rlog(dds_Syl_Vac, blind = FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$treatment, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

######
plotPCA(vsd, intgroup=c("treatment", "Run"))
########
pcaData <- plotPCA(vsd, intgroup=c("treatment", "Run"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=treatment)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

###outlier
W <- res_Syl_Vac$stat
maxCooks <- apply(assays(dds_Syl_Vac)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds_Syl_Vac)
p <- 3
abline(h=qf(.99, p, m - p))
#Top genes : sort by pvalue
resSort <- res_Syl_Vac[order(res_Syl_Vac$pvalue),]
head(resSort)

####
dds_Syl_Vac <- estimateSizeFactors(dds_Syl_Vac)
colSums(counts(dds_Syl_Vac))
dds_Syl_Vac <- dds_Syl_Vac[rowSums(counts(dds_Syl_Vac)) > 0,]
#************************************************************
liver_DESeq_Syl_Vac <- write.table( res_Syl_Vac, "/.../liver_DESeq_Syl_Vac.csv")

hist(res_Syl_Vac$pvalue[res_Syl_Vac$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


DGE.results <- results(dds_Syl_Vac,
                       independentFiltering = TRUE,
                       alpha = 0.05)
table(DGE.results$padj < 0.05)

logcounts_Syl_Vac <- log2( counts(dds_Syl_Vac, normalized=TRUE) + 1 )
pc <- prcomp( t( logcounts_Syl_Vac ) )
plot(hclust(dist(t(logcounts_Syl_Vac))), labels=colData(dds_Syl_Vac)$protocol)
plot(hclust(dist(t(logcounts_Syl_Vac))), labels=colData(dds_Syl_Vac)$treatment)
plot(logcounts_Syl_Vac[,1], logcounts_Syl_Vac[,2], cex=.1)
plotMA(dds_Syl_Vac, ylim=c(-5,5))
plotMA(res_Syl_Vac, ylim=c(-5,5))

###########################################
