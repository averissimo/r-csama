
#
#
# 2:

library(IRanges)
ir <- IRanges(start=c(10, 20, 30), width=5)
ir

flank(ir, 3)
flank(ir, 3, start = FALSE)

class(ir)

getClass(class(ir))

?"flank,Ranges-method" 

#

library(GenomicRanges)
gr <- GRanges(c("chr1", "chr1", "chr2"), ir, strand=c("+", "-", "+"))
gr

flank(gr, 3)

getClass(gr)
class(gr)

?"flank,GenomicRanges-method"

#

showMethods(class="GRanges", where=search())

#

help(package="GenomicRanges")
vignette(package="GenomicRanges")
vignette(package="GenomicRanges", "GenomicRangesHOWTOs")
















