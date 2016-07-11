#
#
#
# Part 1

library(Biostrings)
library(GenomicRanges)

# Documentation of packages
# vignette(package="GenomicRanges")
# vignette("GenomicRangesIntroduction")

# Checkout methods in class
# methods(class="GRanges")

# 1.3

library(Biostrings)                     # Biological sequences
data(phiX174Phage)                      # sample data, see ?phiX174Phage
phiX174Phage

# 1:4 to get only acgt, other positions are for aminoacids
m <- consensusMatrix(phiX174Phage)[1:4,] # nucl. x position counts
polymorphic <- which(colSums(m != 0) > 1)
# check positions where there is a SNP
m[, polymorphic]

methods(class=class(phiX174Phage))      # 'DNAStringSet' methods

#
#
# Exercises

# 1.4
#   1: It is of DNAStringSet class getClass(phiX174Phage)
#      union, consensusMatrix, 

#   2:
vignette(package="Biostrings")
vignette("BiostringsQuickOverview")

#   3: biocViews link is not correct: http://bioconductor.org/packages/release/BiocViews.html#___Software
#      for Biostrings: http://bioconductor.org/packages/release/bioc/html/Biostrings.html

#   4:
library(Biostrings)
data(phiX174Phage)
# count frequency of ACGT in datasets
m <- consensusMatrix(phiX174Phage)[1:4,]
# identify positions of SNPs
polymorphic <- which(colSums(m != 0) > 1)
# show the variations in the SNP positions from each dataset
#  apply substr function (with x = phiX... and subsequent argumets polymorphic, polymorhic)
#  in other words substr gets position in sequence from each database
mapply(substr, polymorphic, polymorphic, MoreArgs=list(x=phiX174Phage))
   
#
#
#
# Part 2

library(GenomicRanges)
library(IRanges)
ir <- IRanges(start=c(10, 20, 30), width=5)
ir
help(package="GenomicRanges")
vignette(package="GenomicRanges")
vignette(package="GenomicRanges", "GenomicRangesHOWTOs")
?GRanges