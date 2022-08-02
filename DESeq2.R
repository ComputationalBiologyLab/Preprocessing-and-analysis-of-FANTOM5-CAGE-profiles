##DESeq2 Script
### Bioconductor and CRAN libraries used
#library(tidyverse)
#library(RColorBrewer)
library(DESeq2)
#library(DEGreport)


##read data and metadata
data = read.csv('F:/MSCs to NSCs/DEGs/NSCs/GeneExpMatrix_F2.csv', 
                header = TRUE, row.names = 1)

metadata = read.csv('F:/MSCs to NSCs/DEGs/NSCs/Metadata.csv', 
                    header = TRUE, row.names = 1)
design <- as.formula(~ sampletype)

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix (countData = data,  colData= metadata,  design = design)

## Run analysis
dds <- DESeq(dds)
#dds <- DESeq(dds, minReplicatesForReplace=Inf)
##plot Dispersion
plotDispEsts(dds)


## DESeq2 results

##These two lines give same output
resultsNames(dds)



##To avoid getting NA values for Pvalue 
#res<- results(dds, name="sampletype_target_vs_background", cooksCutoff=FALSE)
##To avoid getting NA values for Pvalue and Padj
#resnofilt <- results(dds, name="sampletype_target_vs_background", cooksCutoff=FALSE, independentFiltering=FALSE)




##Change the adjusted p value cutoff to be 0.05 and see How many adjusted p-values were less than 0.05?

res05 <- results(dds,name="sampletype_Target_vs_Background", alpha=0.05)
sum(res05$padj < 0.05, na.rm=TRUE)

##Plot Log2fold change
plotMA(res05, ylim=c(-2,2))
write.csv( res05, file="F:/MSCs to NSCs/DEGs/NSCs/Res05_cond2.csv")


##Shrink Log2Fold Change 
res.ape <- lfcShrink(dds=dds, coef=2, type="apeglm")
res.ash <- lfcShrink(dds=dds, coef=2, type="ashr")
res.norm <- lfcShrink(dds=dds, coef=2, type="normal")

par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-2,2)
plotMA(res.norm, xlim=xlim, ylim=ylim, main="normal")
plotMA(res.ape, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(res.ash, xlim=xlim, ylim=ylim, main="ashr")



##HeatMap of CountMatrix
library("pheatmap")
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("group")])
rownames(df) <- colnames(assay(ntd[select,]))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
