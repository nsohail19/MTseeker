library(Sushi) # new dep?
library(MTseeker) 
library(GenomicAlignments) # needs GenomicAlignments::coverage for this step
#MTreads <- readRDS("bulkATAC_MTreads.rds") 

#' plot read coverage across the mitochondrial genome 
#' 
#' @param  mal      an MAlignments (NOT A LIST) 
#' 
#' @import GenomicAlignments
#' @import Sushi 
#' 
#' @export
mtCoveragePlot <- function(mal, bases=NULL) { 

  if (is.null(bases)) {
    data(mtAnno.rCRS)
    bases <- seqlengths(mtAnno)["chrM"]
  }

  # how many reads are on each base?
  asRle <- coverage(as(mal, "GAlignments"))$chrM

  # Sushi expects a data.frame; we can get this from a GAlignments, sort of:
  library(Sushi)
  forSushi <- data.frame(
    chrom=rep("chrM", bases),
    start=seq(1, bases, 1),
    end=seq(1, bases, 1),
    value=as(asRle, "integer")
  )

  plotBedgraph(forSushi, chrom=unique(forSushi$chrom), 
               chromstart=min(forSushi$start), 
               chromend=max(forSushi$end),
               colorbycol=SushiColors(6))
  labelgenome(chrom=unique(forSushi$chrom),
              chromstart=min(forSushi$start), 
              chromend=max(forSushi$end),
              n=floor(bases/1e3), scale="Kb")
  axis(side=2, las=2, tcl=0.2)

}

#bulk1_reads <- MTreads[[1]]

#par(mfrow=c(1,2))
#mtCoveragePlot(bulk1_reads)
#title("Bulk 1")
