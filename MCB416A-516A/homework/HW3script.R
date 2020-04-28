#install DESeq2 at the first time use
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# load DESeq2 (each time to use it)
library(DESeq2)

# the input count table "/Users/izzma/Documents/Counts_Ca.Ctrl_LEfiltered (1).txt" should be in your current working directory (or you can provide the full path to it)
readcounts <- read.table("/Users/izzma/Documents/Counts_Ca.Ctrl_LEfiltered (1).txt", header=TRUE)

# Save $NAME info (col 1) and then remove all columns that don't contain read counts(1)
row.names(readcounts) <- readcounts$NAME 
readcounts <- readcounts[,2:14]
dim(readcounts) 
head(readcounts, n=3) 

# Simplify sample/column names
names(readcounts) <- gsub("_Aligned.sortedByCoord.out.bam", "", names(readcounts)) 
head(readcounts, n=3) 

# make a data frame (to be assigned to colData) with meta-data where row.names should match the individual sample names (i.e., column names in readcounts)
sample_info <- data.frame (row.names=names(readcounts), condition=gsub("_[0-9]+", "", names(readcounts))) 
sample_info 

######
# generate the DESeqDataSet (with required fields: countData, ColData, design)
DESeq.ds <- DESeqDataSetFromMatrix (countData = readcounts, colData = sample_info, design = ~condition) 

# DESeq2's default method to normalize read counts to account for differences in sequencing depths is implemented in estimateSizeFactors()
# calculate the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors (DESeq.ds)

# if you check colData() again, you can see it now contains the sizeFactors
colData (DESeq.ds)

# DESeq2's rlog() function shrinks the variance of low read counts and returns values that are both normalized for sequencing depth and transformed to the log2 scale
DESeq.rlog <- rlog (DESeq.ds, blind = TRUE) 
rlog.counts <- assay(DESeq.rlog) 
head(rlog.counts)

colData(DESeq.ds) 
colData(DESeq.ds)$condition

# set Ctrl as the first-level reference (to be used as the denominator for the fold change calculation) 
colData(DESeq.ds)$condition <- relevel (colData(DESeq.ds)$condition, "Ctrl")
colData(DESeq.ds)$condition  

DESeq.ds <- DESeq (DESeq.ds) # output = DESeqDataSet object
#Note that the input for the DEseq () function are the raw read counts 
#(non-normalized, untransformed), as the function will perform normalizations 
#and transformations under the hood; supplying anything but raw read counts 
#will result in nonsensical results.

#The results() function extracts the base means across samples, moderated log2 fold changes, standard errors, test statistics etc. for every gene.
DGE.results <- results (DESeq.ds, independentFiltering = TRUE , alpha = 0.05) 
summary (DGE.results) 

#check the content of the DGE.results object (like a data frame)
head (DGE.results, n=10)
table (DGE.results$padj < 0.05)
row.names (subset(DGE.results, padj < 0.05))

#write DGE data to output file
write.table(DGE.results, file="DGE.Ca_Ctrl.txt", sep="\t", quote=FALSE, row.names=TRUE)

#make a PCA plot using the plotPCA function 
library (ggplot2)
plotPCA (DESeq.rlog)

# Heatmap
## generating heatmaps using NMF::aheatmap()
#install the NMF package at the first time use
BiocManager::install("NMF")

# load the library with the aheatmap() function
library (NMF) 

# aheatmap() needs a matrix of values: in this example, a matrix of DE genes with the transformed read counts for each replicate
# identify genes with the desired adjusted p-value cut-off
DGEgenes <- row.names(subset(DGE.results, padj < 0.05))

# extract the normalized read counts for DE genes into a matrix
mat_DGEgenes <- rlog.counts[DGEgenes, ]

# plot the normalized read counts of DE genes 
library (NMF) 
aheatmap (mat_DGEgenes,
          Rowv = TRUE , Colv = TRUE,
          distfun="euclidean", hclustfun="median",
          scale="row")

#Find Kegg pathways for upregulated genes
lastq<-subset(DGE.results, log2FoldChange>=1 & padj<0.05)
write.table(lastq, file="genelistHW3.txt", sep="\t", quote=FALSE, row.names=TRUE)

#Filling in the table
table1<-subset(DGE.results, padj<0.05)
write.table(table1, file="table1.txt", sep="\t", quote=FALSE, row.names=TRUE)
table2<-subset(DGE.results, log2FoldChange>=1)
write.table(table2, file="table2.txt", sep="\t", quote=FALSE, row.names=TRUE)
table3<-subset(DGE.results, log2FoldChange<=-1)
write.table(table3, file="table3.txt", sep="\t", quote=FALSE, row.names=TRUE)
