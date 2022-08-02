x <- read.delim("/Users/cesc/Desktop/Asmaa/MSCs_Target/human_GRCh38_mart.txt", stringsAsFactors = FALSE) 
e <- GRanges(paste0("chr", x$Chromosome.scaffold.name),
             IRanges(x$Exon.region.start..bp., x$Exon.region.end..bp.),
             ifelse(x$Strand + 1, "+", "-"))
e$gene_name <- Rle(x$Gene.name)
e$transcript_type <- Rle(x$Gene.type) 
e$type <- "exon" 
e$type <- Rle(e$type)


##Exon Object
e <- GRanges(paste0("chr", x$Chromosome.scaffold.name), 
             IRanges(x$Exon.region.start..bp., x$Exon.region.end..bp.),
             ifelse(x$Strand + 1, "+", "-")) 
e$gene_name <- Rle(x$Gene.name) 
e$transcript_type <- Rle(x$Gene.type) 
e$type <- "exon" 
e$type <- Rle(e$type) 
e <- sort(unique(e))


##Gene Object
g <- GRanges( paste0("chr", x$Chromosome.scaffold.name)
                 , IRanges(x$Gene.start..bp., x$Gene.end..bp.)
                 , ifelse( x$Strand + 1, "+", "-"))
g$gene_name <- Rle(x$Gene.name) 
g$transcript_type <- Rle(x$Gene.type) 
g$type <- "gene" 
g$type <- Rle(g$type) 
g <- sort(unique(g))

##Transcript Object

t <- GRanges( paste0("chr", x$Chromosome.scaffold.name)
                 , IRanges(x$Transcript.start..bp., x$Transcript.end..bp.)
                 , ifelse( x$Strand + 1, "+", "-"))



t$gene_name <- Rle(x$Gene.name) 
t$transcript_type <- Rle(x$Gene.type) 
t$type <- "transcript" 
t$type <- Rle(t$type) 
t <- sort(unique(t))



gff <- sort(c(g, t, e))
seqlevels(gff) <- seqlevelsInUse(gff)

save(gff , file = "human_ann.RData", compress = "xz")
