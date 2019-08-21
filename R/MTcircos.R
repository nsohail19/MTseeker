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
                     incol=NULL, anno=NULL, how=c("matrix", "AA"), ...) {
  circos.clear() 

  if (length(how) > 1) {
    how <- "matrix"
  }

  if (how == "AA" && !("typeMut" %in% names(mcols(variants[[1]])))) {
    message("Must run decomposeAndCalcConsequences() before plottting AA changes")
    stop()
  }
  
  anno <- initMTcircos(variants)
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  
  if (!is.null(variants) && length(variants) != 0) {
    message("Splitting variants by strand...")
    stranded <- byStrand(variants, anno)
    message("Replacing `outside` with heavy-strand variants...")
    outside <- stranded$heavy
    message("Replacing `inside` with light-strand variants...")
    inside <- stranded$light
  } 

  # Outside track
  if (!is.null(outside) && length(outside) != 0) {
    
    # Color code according to AA changes
    if (how == "AA") {
      message("Working on it, AA")
    }
    
    else {
      # Color coding for the variants
      # del=blue, SNV=black, ins=red
      bed1 <- .makeColoredMatrix(outside)
      vafBed1 <- .vafMatrix(bed1, outside)

      #col_fun = colorRamp2(c(0, 1, 2, 3), c("white", "red", "black", "blue"))
      circos.genomicHeatmap(bed1, outcol, line_col=.colorCode(bed1$chr), col = vafBed1[,4:ncol(vafBed1)],
                            track.margin=c(0,0), side="outside", border=NA,
                            line_lwd=2) 
      
    }
    
  }
  else {
    circos.track(track.height=0.15, ylim=c(0,1), bg.border=NA)
  }
  
  # main track, gene names and such
  res <- genesMTcircos(variants, anno, legends=T)
  
  # Inside track
  if (!is.null(inside) && length(inside) != 0) {
    
    # Color code according to AA changes
    if (how == "AA") {
      message("Working on it, AA")
    }
    
    else {
      # del=blue, SNV=black, ins=red

      bed2 <- .makeColoredMatrix(inside)
      vafBed2 <- .vafMatrix(bed2, inside)

      #col_fun = colorRamp2(c(0, 1, 2, 3), c("white", "red", "black", "blue"))
      circos.genomicHeatmap(bed2, outcol, line_col=.colorCode(bed2$chr), col = vafBed2[,4:ncol(vafBed2)],
                            track.margin=c(0,0), side="inside", border=NA,
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
  if (is(mvr, "MVRangesList")) {
    allNames <- lapply(mvr, names)
    allNames <- unique(unlist(unname(allNames)))
    numSamples <- length(mvr)
  } else {
    allNames <- unique(names(mvr))
    numSamples <- 1
  }

  rowNam <- c("chr", "start", "end")

  if (numSamples == 1) rowNam <- append(rowNam, unique(as.character(sampleNames(mvr))))
  else rowNam <- append(rowNam, names(mvr))
  
  m <- matrix(0, ncol = length(rowNam), nrow = length(allNames))
  typeDF <- data.frame(m)

  names(typeDF) <- rowNam
  rownames(typeDF) <- allNames

  # Figure out which of the variants in a sample is SNV, ins, del
  # Assign color according to the unique values set for each type of variant
  for (i in 1:numSamples) {
    
    # MVRangesList
    if (! (numSamples == 1) ) {
      rowInd <- which(allNames %in% names(mvr[[i]]))
      rowOverlapNames <- allNames[rowInd]
    }
    # MVRanges
    else {
      rowOverlapNames <- allNames
    }

    snvs <- rowOverlapNames[grep(">", rowOverlapNames)]
    ins <- rowOverlapNames[grep("ins", rowOverlapNames)]
    dels <- rowOverlapNames[grep("del", rowOverlapNames)]
    
    # del=blue, SNV=black, ins=red
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

  # Get the start and end position for each variants
  for (j in 1:numSamples) {

    hits <- which(typeDF[,j + 3] != 0)
    varNames <- allNames[hits]
    
    vars <- mvr[[j]][varNames]
    varStart <- start(vars)
    varEnd <- end(vars)
    
    typeDF$start[hits] <- varStart
    typeDF$end[hits] <- varEnd
  }

  # Get the names of the genes each variant is found in
  anno <- suppressMessages(getAnnotations(mvr))
  ov <- findOverlaps(IRanges(typeDF$start, typeDF$end), ranges(anno))

  # If there are overlapping genes
  # Only list the first gene it overlaps in
  firstOv <- ov[match(unique(queryHits(ov)), queryHits(ov)), ]
  
  typeDF$chr <- names(anno)[subjectHits(firstOv)]

  return(typeDF)
}

.vafMatrix <- function(bed, mvr) {

  vafs <- bed
  vafs[,4:ncol(vafs)] <- 0

  # Get the VAF for each variants
  for (j in 1:(ncol(bed) - 3)) {
    
    hits <- which(bed[,j + 3] != 0)
    varNames <- row.names(bed)[hits]
    
    vars <- mvr[[j]][varNames]
    varVAF <- vars$VAF
    
    vafs[,j + 3][hits] <- varVAF
  }

  # Transparancy goes from a scale of 0-255
  #copy[,4:ncol(vafs)] * 255
  #return(vafs)
  
  #vafBed1 <- .vafMatrix(typeDF, mvr)
  #adjustcolor(typeDF[,4][1], alpha.f=vafBed1[,4][1])

  # Nonzero elements
  
  nonzero <- which(bed[,4:ncol(bed)] !=0, arr.ind=T)
  for (k in 1:nrow(nonzero)) {
    
    # 1 is the row
    # 2 is the column
    rows <- nonzero[k,][1]
    cols <- nonzero[k,][2]
    
    #col_fun = colorRamp2(c(0, 1, 2, 3), c("white", "red", "black", "blue"))
    if (bed[rows,][cols + 3] == 1) bed[rows,][cols + 3] <- "red"
    else if (bed[rows,][cols + 3] == 2) bed[rows,][cols + 3] <- "black"
    else if (bed[rows,][cols + 3] == 3) bed[rows,][cols + 3] <- "blue"

    bed[rows,][cols + 3] <- adjustcolor(col = bed[rows,][cols + 3], 
                                        alpha.f = vafs[rows,][cols + 3])
  }
   
  return(bed)
}
