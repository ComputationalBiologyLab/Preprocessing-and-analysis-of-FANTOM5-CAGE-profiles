library(CAGEr)

##Path to raw data

pathsToInputFiles <- list.files("/Users/cesc/Desktop/Asmaa/MSCs_Target",
                                pattern = "*.ctss.bed",      full.names = TRUE)


##Create ce object of CAGEr class

ce <- CAGEexp( genomeName     = "BSgenome.Hsapiens.UCSC.hg38"
               , inputFiles     = pathsToInputFiles
               , inputFilesType = "ctss"
               , sampleLabels   = sub( ".ctss.bed", "", basename(pathsToInputFiles)))

##To read the data 
getCTSS(ce)


CTSStagCountSE(ce)

CTSScoordinatesGR(ce)
CTSStagCountDF(ce)
head(CTSScoordinates(ce))
head(CTSStagCountDf(ce))
head(CTSStagCount(ce))
sampleLabels(ce)


##Quality Control
#(1) Annotation

library(SummarizedExperiment)
load("/Users/cesc/Desktop/Asmaa/MSCs_Target/human_ann.RData")
annotateCTSS(ce, gff)
colData(ce)[,c("librarySizes", "promoter", "exon", "intron", "unknown")]
plotAnnot(ce, "counts")


#Add a gene expression table in the GeneExpSE experiment slot of an annotated CAGEexp object. 
CTSStoGenes(ce)
all( librarySizes(ce) -
                   colSums(SummarizedExperiment::assay(GeneExpSE(ce))) ==
                   ce$unannotated)

##Summary about the ce after conversion to gene symbol
GeneExpSE(ce)
#GeneExp Matrix
GeneExpMatrix <- assay(GeneExpSE(ce))

