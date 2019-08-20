#' plot a canonical human (or, in principle, any) mitochondrial genome 
#'
#' The default font sizes, orientations, etc. are optimized for a "cold" start;
#' if you want to fiddle with the details, crack open the code and modify it...
#' or alternatively, add sectors/dendrograms inside of this "framed" version.
#'
#' 
#' @param variants  optional MVRanges or MVRangesList to split by strand & plot
#' @param outside   optional MVRanges or MVRangesList to plot outside the circle
#' @param inside    optional MVRanges or MVRangesList to plot inside the circle
#' @param outcol    optional color assignment function or matrix for outside
#' @param incol     optional color assignment function or matrix for inside
#' @param anno      a GRanges (optional, defaults to mtAnno.rCRS if none given)
#' @param how       optional specification for how to plot multiple samples
#' @param ...       other arguments to pass on to called functions
#' 
#' @return          invisibly, a list: `anno` (data.frame) + `pfun` (panel.fun)
#'
#' @import VariantTools
#' @import rtracklayer
#' @import circlize
#' @import viridis
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @examples 
#' library(MTseekerData)
#' data(RONKSvariants) 
#' MTcircos(RONKSvariants)
#' # same as plot(RONKSvariants)
#' title("Renal oncocytomas and normal kidney samples")
#' 
#' @export 
MTcircos <- function(variants=NULL, AA=FALSE, outside=NULL, inside=NULL, outcol=NULL, 
                     incol=NULL, anno=NULL, how=c("matrix","VAF", "AA"), ...) {
  circos.clear() 
  
  if (length(how) > 1) {
    how <- "matrix"
  }
  
  if (how == "AA") {
    message("Must run decomposeAndCalcConsequences() before plottting AA changes")
    stop()
  }
  
  anno <- initMTcircos(variants)
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  
  if (!is.null(variants)) {
    message("Splitting variants by strand...")
    stranded <- byStrand(variants, anno)
    message("Replacing `outside` with heavy-strand variants...")
    outside <- stranded$heavy
    message("Replacing `inside` with light-strand variants...")
    inside <- stranded$light
  } 
  
  # Outside track
  if (!is.null(outside)) {
    
    # Color code according to AA changes
    if (how == "AA") {
      message("Working on it, AA")
    }
    
    # Transparency according to VAF
    else if (how == "VAF") {
      message("Working on it, VAF")
    }
    
    else {
      # Color coding for the variants
      # del=blue, SNV=black, ins=red
      bed1 <- .makeColoredMatrix(outside)
      col_fun = colorRamp2(c(0, 1, 2, 3), c("white", "red", "black", "blue"))
      circos.genomicHeatmap(bed1, outcol, line_col=.colorCode(bed1$chr), col = col_fun,
                            track.margin=c(0,0), side="outside", border=NA,
                            line_lwd=2) 
      
    }
    
  }
  else {
    circos.track(track.height=0.15, ylim=c(0,1), bg.border=NA)
  }
  
  # main track, gene names and such
  res <- genesMTcircos(variants, anno)
  
  # Inside track
  if (!is.null(inside)) {
    
    # Color code according to AA changes
    if (how == "AA") {
      message("Working on it, AA")
    }
    
    # Transparency according to VAF
    else if (how == "VAF") {
      message("Working on it, VAF")
    }
    
    else {
      # Color coding for the variants
      # del=blue, SNV=black, ins=red
      bed2 <- .makeColoredMatrix(inside)
      col_fun = colorRamp2(c(0, 1, 2, 3), c("white", "red", "black", "blue"))
      circos.genomicHeatmap(bed2, outcol, line_col=.colorCode(bed1$chr), col = col_fun,
                            track.margin=c(0,0), side="outside", border=NA,
                            line_lwd=2) 
      
    }
    
  }
  else {
    circos.track(track.height=0.15, ylim=c(0,1), bg.border=NA)
  }
  
  # Color code according to AA changes
  if (how == "AA") {
    message("Working on it, AA")
  }
  
  # Transparency according to VAF
  else if (how == "VAF") {
    message("Working on it, VAF")
  }
  
  else {
    legend("bottomleft", title="Type of Variant",
           legend=c("Insertion", "Deletion", "SNV"), col=c("red", "blue", "black"), pch=15, cex=0.8)
  }
 
  
  invisible(res)
  
}


# helper fn
.makeBed <- function(x) {
  bed <- switch(class(x),
                "MVRanges"=.mvrToBed(x),
                "MVRangesList"=.mvrlToBed(x),
                "GRanges"=.grToBed(x))
  names(bed)[1] <- "chr"
  return(bed)
}

# helper fn
.mvrToBed <- function(mvr) { 
  
  message("This will take a moment")

  # Iterate through each variant and call locateVariant
  newMvr <- locateMTvariants(mvr[1])
  for (i in 2:length(mvr)) {
    newVar <- locateVariants(mvr[i])
    newMvr <- append(newMvr, newVar)
  }
  
  bed <- as.data.frame(newMvr)[, c("gene", "start", "end")]
  
  bed$value <- mvr$VAF
  bed <- subset(bed, !is.na(bed[,1]))
  return(bed)
}

# helper fn
.mvrlToBed <- function(mvrl) { 

  gr <- granges(mvrl, filterLowQual=FALSE)
  bed <- as.data.frame(gr)
  bed <- bed[, c(6, 2, 3, 8:(ncol(bed)))]
  return(bed)
}

# helper fn
.grToBed <- function(gr) { 
  bed <- as.data.frame(gr)[, c("gene", "start", "end")]
  bed$value <- 1
  return(bed)
}


.makeColoredMatrix <- function(mvr) {
  
  # Making a colored matrix
  allNames <- lapply(mvr, names)
  allNames <- unique(unlist(unname(allNames)))
  
  rowNam <- c("chr", "start", "end")
  rowNam <- append(rowNam, names(mvr))
  
  m <- matrix(0, ncol = length(rowNam), nrow = length(allNames))
  typeDF <- data.frame(m)
  names(typeDF) <- rowNam
  rownames(typeDF) <- allNames
  
  for (i in 1:length(mvr)) {
    
    rowInd <- which(allNames %in% names(mvr[[i]]))
    rowOverlapNames <- allNames[rowInd]
    
    snvs <- rowOverlapNames[grep(">", rowOverlapNames)]
    ins <- rowOverlapNames[grep("ins", rowOverlapNames)]
    dels <- rowOverlapNames[grep("del", rowOverlapNames)]
    
    if (length(ins)) {
      typeDF[ins,][,i + 3] <- 1
    }
    
    if (length(snvs)) {
      typeDF[snvs,][,i + 3] <- 2
    }
    
    if (length(dels)) {
      typeDF[dels,][,i + 3] <- 3
    } 
  }
  
  split <- gsub("[^0-9\\_]", "", allNames) 
  pos <- as.numeric(sub("_[^_]+$", "", split))
  
  typeDF$start <- pos
  typeDF$end <- pos
  
  anno <- suppressMessages(getAnnotations(testPu_anno))
  ov <- findOverlaps(IRanges(typeDF$start, typeDF$end), ranges(anno))
  
  typeDF$chr <- names(anno)[subjectHits(ov)]
  
  return(typeDF)
}