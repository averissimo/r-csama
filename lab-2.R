#' ---
#' title: "R Markdown first test"
#' output: 
#'   github_document
#' date: "It's today!! duhhhh"
#' ---

#' load library with the data
library("airway")

#' base dir with all the data files (csv)
dir <- system.file("extdata", package="airway", mustWork=TRUE)
dir

list.files(dir)

#' table with overview of the existing files
csvfile <- file.path(dir,"sample_table.csv")
sampleTable <- read.csv(csvfile,row.names=1)
sampleTable

#' set a vector with existing file paths
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))

#' # file.exists(filenames)

#' interface to BAM files
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

seqinfo(bamfiles[8])

#' loading pre-built transcript databse (GTF files)
library("GenomicFeatures")
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")

#' exons grouped by gene
ebg <- exonsBy(txdb, by="gene")
ebg

#' register more cores
library("BiocParallel")
register(MulticoreParam(3))

#' actual count 
library("GenomicAlignments")
se <- summarizeOverlaps(features=ebg,       # exons by gene
                        reads=bamfiles,     # bamfiles interface
                        mode="Union",       # map if it partial aligns with at most 1 gene
                        singleEnd=FALSE,    # paired end reads
                        ignore.strand=TRUE, # experiment was not strand specific
                        fragments=TRUE )    # count fragments if unpaired

#' 'actual view of counts in summarized experiment object (se)
head(assay(se))

dim(se) # 20 genes and 8 samples

rowRanges(se)

rowRanges(se)[[4]]

#' 'display metadata (shows origin of analysis)
str(metadata(rowRanges(se)))

colData(se)

colnames(se) == paste0(sampleTable$Run, '_subset.bam')
colnames(se) == names(bamfiles)

colData(se) <- DataFrame(sampleTable)
colData(se)

#' '
#' '
#' 'Alternative to process above for counting
library("Rsubread")
fc <- featureCounts(files=filenames, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE)

colnames(fc$counts) <- sampleTable$Run

#' 'test if the counts are identical
ixs_dif <- which(!sapply(1:nrow(se), function(ix){ return(all(assay(se)[ix,])) }))
assay(se)[rownames(fc$counts)[ixs_dif],]
fc$counts[rownames(fc$counts)[ixs_dif],]
if (length(ixs_dif) > 0)
  print(paste0('Counts are not the same for: ', paste(rownames(fc$counts)[ixs_dif], collapse = ' ')))
#' 'they are not!

#' 'end of alternative
#' '

colData(se)

#' 're-order level to make untrt first
se$dex
se$cell 

se$dex <- relevel(se$dex, "untrt")
se$dex

#' '
#' '
#' '
data("airway")
se <- airway # use the publicly available summarized experiment from airway package

#' 'pre-filter with only counts higher than 5
se <- se[ rowSums(assay(se)) >= 5, ]

#' 're-order level
se$dex <- relevel(se$dex, "untrt")
se$dex

colData(se)

#' 'Build a DESeq2 specific summarized experiments
#' called DESeqDataSet
library(DESeq2)
dds <- DESeqDataSet(se, design = ~ cell + dex)

#' 'to build a DESeqDataSet from count matrix one would do
countdata <- assay(se)
head(countdata, 3)
#' '
coldata <- colData(se)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)

#' '
#' '
#' 'now for edgeR

library("edgeR")
genetable <- data.frame(gene.id=rownames(se))
y <- DGEList(counts=countdata, 
             samples=coldata, 
             genes=genetable)
names(y)

#' '
#' 'PCA
#' '
#' 'normalize the count to allow for PCA analysis
vsd <- vst(dds)

head(assay(vsd),3)
head(assay(dds),3)

#' 'plot PCA using default options
plotPCA(vsd, "dex")

#' 'improve on the plot by using ggplot
data <- plotPCA(vsd, intgroup = c( "dex", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
library("ggplot2")
ggplot(data, aes(PC1, PC2, color=dex, shape=cell)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

#' '
#' '
#' 'MDS (similar analysis as PCA)
#' multidimensional scaling
#' '
y <- calcNormFactors(y)
plotMDS(y, top = 1000, labels = NULL, col = as.numeric(y$samples$dex), 
        pch = as.numeric(y$samples$cell), cex = 2)

#' using ggplot2 and manual distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(vsd)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)

#' using poisson distances
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
ggplot(mdsPois, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)
  
#' 
#' 
#' 
#' Differential expression analysis
#' 

dds <- DESeq(dds)
#' for p-values < .1
res <- results(dds)
mcols(res, use.names=TRUE)
summary(res)
table(res$padj < .1)

#' for p-values < .05
res.05 <- results(dds, alpha=.05)
summary(res.05)

#' raising the log2 fold change threshold
#' p-values < .1 and difference is higher than 2 and lower than 1
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < .1)
#' sampleX <- 10^(-1:-20)
#' plot(1:length(sampleX),sapply(sampleX, function(p) {as.vector(table(resLFC1$padj < p)[2])}),
#'     xlab = 'p-values', ylab = 'Count of p-values'); grid()

#' design matrix
#'  simple with just cell line and treated/untreated

y <- DGEList(counts=countdata, 
             samples=coldata, 
             genes=genetable)

design <- model.matrix(~ cell + dex, y$samples)
colnames(design)
#' calculate normalization factors
y <- calcNormFactors(y)
#' calculate dispersion
y <- estimateDisp(y, design)

#' calculate GLM
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=ncol(design))
tt <- topTags(lrt, n=nrow(y), p.value=.1)
tt10 <- topTags(lrt, n=10) # just the top 10 by default
tt10

head(lrt$fitted.values)
head(y$counts)
head(deviance(lrt))

#' compare between two software overlap (DESeq2 and edgeR)
tt.all <- topTags(lrt, n=nrow(y), sort.by="none")
table(DESeq2=res$padj < 0.1, edgeR=tt.all$table$FDR < 0.1)

#' compare for fold-change threshold
treatres <- glmTreat(fit, coef = ncol(design), lfc = 1)
tt.treat <- topTags(treatres, n = nrow(y), sort.by = "none")
table(DESeq2 = resLFC1$padj < 0.1, edgeR = tt.treat$table$FDR < 0.1)

#' rank them
common <- !is.na(res$padj)
plot(rank(res$padj[common]), 
     rank(tt.all$table$FDR[common]), cex=.1,
     xlab="DESeq2", ylab="edgeR"); grid()

#' 
#' 
#' 
#' Plotting the results
#' 

topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, topGene, "dex")

#' MA plot for DESeq2
DESeq2::plotMA(res, ylim=c(-5,5))
#' MA plot for edgeR
plotSmear(lrt, de.tags=tt$table$gene.id)

#' 
#' Heatmap of 30 most significant genes
library("pheatmap")
mat <- assay(vsd)[ head(order(res$padj),30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)

#' 
#' 
#' 
#' Annotate genes with gene name
#' 
library("AnnotationDbi")
library("Homo.sapiens")

columns(Homo.sapiens)

#' using column ensembl to map genes
#'  attention: it will show a message that should be ignored!
res$symbol <- mapIds(Homo.sapiens,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
y$genes$symbol <- res$symbol

res$entrez <- mapIds(Homo.sapiens,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

y$genes$entrez <- res$entrez

res$genename <- mapIds(Homo.sapiens,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")

y$genes$genename <- res$genename

resOrdered <- res[order(res$padj),]
head(resOrdered)

#' export results
resOrderedDF <- as.data.frame(resOrdered)[seq_len(100),]
write.csv(resOrderedDF, file="results.csv")

#' prettify it with html results!

library("Glimma")
glMDPlot(lrt, 
         counts=y$counts, 
         anno=y$genes, 
         groups=y$samples$dex, 
         samples=colnames(y),
         status=tt.all$table$FDR < 0.1,
         id.column="gene.id")

#' DESeq2 results
res.df <- as.data.frame(res)
res.df$log10MeanNormCount <- log10(res.df$baseMean)
idx <- rowSums(counts(dds)) > 0
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1
glMDPlot(res.df,
         xval="log10MeanNormCount",
         yval="log2FoldChange",
         counts=counts(dds)[idx,],
         anno=data.frame(GeneID=rownames(dds)[idx]),
         groups=dds$dex,
         samples=colnames(dds),
         status=res.df$padj < 0.1,
         display.columns=c("symbol", "entrez"))

#' continue at: Gene set overlap analysis
