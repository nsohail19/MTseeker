locateMTvariants <- function(query) {

  # Checks if the function has been run before
  if ("gene" %in% names(mcols(query)) &
      "region" %in% names(mcols(query)) &
      "localEnd" %in% names(mcols(query)) & 
      "localStart" %in% names(mcols(query)) &
      "startCodon " %in% names(mcols(query)) &
      "endCodon" %in% names(mcols(query))) {
    return(query) # done 
  }
  
  # No input was given
  if (length(query) == 0) return(NULL)
  
  #run serially
  if (is(query, "MVRangesList")) {
    mvr <- MVRangesList(lapply(query, locateVariants))
    return(mvr)
  }
  
  # For now only support rCRS
  stopifnot(genome(query) == "rCRS")
  data("mtAnno.rCRS", package="MTseeker")
  metadata(query)$annotation <- mtAnno
  
  # Localized genic coordinates
  # Ben said he was interested in looking at tRNA regions 
  #anno <- subset(mtAnno, region %in% c("tRNA", "coding"))
  
  # decomposeAndCalc throws me an error when I try to use tRNA
  anno <- subset(mtAnno, region == "coding")
  ol <- findOverlaps(query, anno, ignore.strand=TRUE)
  
  if (length(ol) == 0) {
    message("No overlapping genes in coding regions for variant: ", names(query))
    return(query[0])
  }
  
  if (length(subjectHits(ol)) > 2) {
    warning("More than 2 overlapping genes in locateVariants")
    stop()
  }
  
  query$gene <- NA_character_
  query$overlapGene <- NA_character_

  # If there are multiple genes the variant is located within
  # Create multiple copies of the variant to handle each gene that it is located within
  if (length(subjectHits(ol)) > 1) {

    query <- rep(query, length(subjectHits(ol)))
    
    # Assume that there are only 2 overlapping genes
    query[1]$gene <- names(anno)[subjectHits(ol)][1]
    query[2]$gene <- names(anno)[subjectHits(ol)][2]
    query$overlapGene <- paste(names(anno)[subjectHits(ol)], collapse = ",")
  } 
  
  # Otherwise there is no overlap
  else {
    query$gene <- names(anno)[subjectHits(ol)]
  }
  
  # Initialize
  query$region <- NA_character_
  query$localStart <- NA_integer_
  query$localEnd <- NA_integer_
  query$startCodon <- NA_integer_
  query$endCodon <- NA_integer_
  
  for (i in 1:length(query)) {
    
    subHit <- subjectHits(ol)[i]

    # FOR SOME REASON MT-ND6 DOES NOT CREATE THE TYPICAL 
    if (names(anno[subHit]) == "MT-ND6") {
      
      #browser()
      #MT_CODE <- getGeneticCode("SGC1")
      #gr <- genes(query)
      #refDNA <- getSeq(rCRSeq, gr["MT-ND6"])
      #refAA <- suppressWarnings(translate(refDNA, MT_CODE))

      #offset <- suppressWarnings(translate(extractAt(unlist(refDNA), IRanges(4, width(refAA))), MT_CODE))
    }

    # Does the variant fall into a coding, tRNA, rRNA, or D-loop region
    query[i]$region <- anno[subHit]$region
    
    # Local start and end (relative to gene of interest)
    # 1-indexed
    query[i]$localStart <- start(query[i]) - start(anno[subHit]) + 1
    query[i]$localEnd <- end(query[i]) - start(anno[subHit]) + 1
    
    if (genome(query) == "rCRS" && names(anno[subHit]) == "MT-ND6") {
      
      query[i]$localStart <- start(query[i]) - start(anno[subHit]) + 1 - 2
      query[i]$localEnd <- end(query[i]) - start(anno[subHit]) + 1 - 2
    }
    
    # Local codon starts and ends
    # Taking a conservative approach
    query[i]$startCodon <- floor(query[i]$localStart / 3)
    query[i]$endCodon <- ceiling(query[i]$localEnd / 3)

    # Everything is 1-indexed (if variant is in the first codon)
    if (query[i]$startCodon == 0) query[i]$startCodon <- 1
    if (query[i]$endCodon == 0) query[i]$endCodon <- 1
  }

  return(query)
}
