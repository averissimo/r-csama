# load library with the data
library("airway")

# base dir with all the data files (csv)
dir <- system.file("extdata", package="airway", mustWork=TRUE)
dir

list.files(dir)

# table with overview of the existing files
csvfile <- file.path(dir,"sample_table.csv")
sampleTable <- read.csv(csvfile,row.names=1)
sampleTable

# set a vector with existing file paths
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))

# file.exists(filenames)

# interface to BAM files
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

seqinfo(bamfiles[8])

# loading pre-built transcript databse (GTF files)
library("GenomicFeatures")
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")

# exons grouped by gene
ebg <- exonsBy(txdb, by="gene")
ebg

# register more cores
library("BiocParallel")
register(MulticoreParam(3))

# actual count 
library("GenomicAlignments")
se <- summarizeOverlaps(features=ebg,       # exons by gene
                        reads=bamfiles,     # bamfiles interface
                        mode="Union",       # map if it partial aligns with at most 1 gene
                        singleEnd=FALSE,    # paired end reads
                        ignore.strand=TRUE, # experiment was not strand specific
                        fragments=TRUE )    # count fragments if unpaired

# actual view of counts in summarized experiment object (se)
head(assay(se))

dim(se) # 20 genes and 8 samples

rowRanges(se)

rowRanges(se)[[4]]

# display metadata (shows origin of analysis)
str(metadata(rowRanges(se)))

colData(se)

colnames(se) == paste0(sampleTable$Run, '_subset.bam')
colnames(se) == names(bamfiles)

colData(se) <- DataFrame(sampleTable)
colData(se)

#
#
# Alternative to process above for counting
library("Rsubread")
fc <- featureCounts(files=filenames, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

colnames(fc$counts) <- sampleTable$Run

# test if the counts are identical
ixs_dif <- which(!sapply(1:nrow(se), function(ix){ return(all(assay(se)[ix,])) }))
assay(se)[rownames(fc$counts)[ixs_dif],]
fc$counts[rownames(fc$counts)[ixs_dif],]
if (length(ixs_dif) > 0)
  print(paste0('Counts are not the same for: ', paste(rownames(fc$counts)[ixs_dif], collapse = ' ')))
# they are not!

# end of alternative
#

colData(se)

# re-order level to make untrt first
se$dex
se$cell 

se$dex <- relevel(se$dex, "untrt")
se$dex

#
#
#
data("airway")
se <- airway # use the publicly available summarized experiment from airway package

# pre-filter with only counts higher than 5
se <- se[ rowSums(assay(se)) >= 5, ]

# re-order level
se$dex <- relevel(se$dex, "untrt")
se$dex

colData(se)

# Build a DESeq2 specific summarized experiments
#  called DESeqDataSet
library(DESeq2)
dds <- DESeqDataSet(se, design = ~ cell + dex)

# to build a DESeqDataSet from count matrix one would do
countdata <- assay(se)
head(countdata, 3)
#
coldata <- colData(se)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)

#
#
# now for edgeR

library("edgeR")
genetable <- data.frame(gene.id=rownames(se))
y <- DGEList(counts=countdata, 
             samples=coldata, 
             genes=genetable)
names(y)

#
# PCA
#
# normalize the count to allow for PCA analysis
vsd <- vst(dds)

head(assay(vsd),3)
head(assay(dds),3)

# plot PCA using default options
plotPCA(vsd, "dex")

# improve on the plot by using ggplot
data <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

#
#
# MDS (similar analysis as PCA)
#  multidimensional scaling
#
y <- calcNormFactors(y)
plotMDS(y, top = 1000, labels = NULL, col = as.numeric(y$samples$dex), 
        pch = as.numeric(y$samples$cell), cex = 2)

# using ggplot2 and manual distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)

# using poisson distances
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
ggplot(mdsPois, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)
  
#
#
#
# Differential expression analysis
# 










