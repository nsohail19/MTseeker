#' Annotate AA changes in MT variants
#'
#' @name getProteinImpact
#'
#' @param mvr    An MVRangesList or MVRanges object
#' @param parallel    Whether to run things in parallel
#' @param ...    Other arguments to pass to injectMTVariants
#'
#' @return    Annotated variants
#'
#' @examples
#' 
#' 
#' 
#' 
getProteinImpact <- function(mvr, parallel=FALSE, cores=1) {

  stopifnot(is(mvr, "MVRanges"))
  
  message("Processing consequences for ", sampleNames(mvr)@values)
  
  # Right now can only handle variants in coding regions
  mvr <- .getCodingRegions(mvr)
  
  # Return empty MVRanges if there are none in coding regions
  if (length(mvr) == 0) return(mvr)
  
  # Output that we want to have by the end of this
  mvr$AAchange <- NA_character_
  mvr$impacted.gene <- NA_character_
  mvr$overlapGene <- NA_character_
  mvr$apogee <- TRUE
  
  for (i in 1:length(mvr)) {

    # Get the url and scrape the information from mitimpact
    pos <- paste0(start(mvr[i]),  "-", end(mvr[i]))
    url <- paste0(paste("http://mitimpact.css-mendel.it", "api", "v2.0", 
                        "genomic_position/", sep="/"), pos)
    res <- as.data.frame(read_json(url, simplifyVector=TRUE)$variants, stringsAsFactors=F)
    
    # If mitimpact gives a valid result
    if (nrow(res) > 0 ) {
      
      # Rename the variant and paste together a simpler version of the consequences
      res$genomic <- with(res, paste0("m.", Start, Ref, ">", Alt))
      res$consequence <- with(res, paste0(AA_ref,AA_position,AA_alt))
      
      # The information we are interested in
      columns <- c("genomic","protein","Start", "Ref","Alt","Codon_substitution","consequence",
                   "Mitomap_Phenotype","Mitomap_Status","Gene_symbol",
                   "OXPHOS_complex","Consequence","APOGEE_boost_consensus")
      columnsPresent <- intersect(columns, names(res))
      mitOut <- res[, columnsPresent]
      
      # Sometimes there are multiple variants that can be found at one position
      # Find the one that aligns with the variant we are looking at
      index <- which(mitOut$Alt == alt(mvr[i]))
      
      # If none of the results from apogee have the same alt and ref sequence
      # Go back to the old way of doing things
      if (length(index) == 0) {
        mvr[i] <- decomposeAndCalcConsequences(mvr[i])
        mvr[i]$apogee <- FALSE
      } 
      
      # If apogee gives a result that we want
      else {
        
        if (mitOut$Ref[index] != ref(mvr[i])) {
          message("Reference disparity between APOGEE and pileup")
        }
        
        out <- mitOut[index,]
        mvr[i]$AAchange <- out$consequence
        mvr[i]$impacted.gene <- out$Gene_symbol
      }
    }
    
    # If mitimpact does not give a valid result, go back to the old way
    else {
      mvr[i] <- decomposeAndCalcConsequences(mvr[i])
      mvr[i]$apogee <- FALSE
    }
    
  } #for
 
  # When using decomp, if there are no AA changes, it will set it to ""
  # Setting it to NA if someone wants to remove synonymous mutation
  mvr[which(mvr$AAchange == "")]$AAchange <- NA_character_
  mvr[which(!mvr$apogee)]$impacted.gene <- paste0("MT-", mvr[which(!mvr$apogee)]$impacted.gene) 
  
  mvr <- .typeOfMutation(mvr)
  
  return(mvr) 
}

# helper function to subset ranges to just coding space
.getCodingRegions <- function(mvr) {
  
  # rCRS only, for the time being 
  stopifnot(unique(genome(mvr)) == "rCRS")
  
  # get mtGenes if needed 
  gr <- genes(mvr)
  stopifnot(unique(genome(gr)) == "rCRS")
  
  mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within"))
  
  # subset the variants to those that overlap the target GRanges and are canon
  if (length(mvr)) {
    
    #drop anything that has an N base.. this also looks like a weird bug?
    mvr <- mvr[!grepl("N", mvr@alt),]
    
    #check again whether we've now cleared out all the variants
    #return an empty ranges if we have
    if (length(mvr) == 0) {
      mvr <- MVRanges(subsetByOverlaps(mvr, gr, type="within"))
    }
  } 
  return(mvr)
}

# helper fn to determine type of mutation
.typeOfMutation <- function(mvr) {
  
  mvr$typeMut <- NA_character_

  # Determine the overaching consequences (type of mutation)
  # SNP
  snp <- grep(">", names(mvr))
  if (length(snp) != 0) {
    
    # Synonymous : no AA change overall
    syn <- which(is.na(mvr[snp]$AAchange))
    mvr[snp]$typeMut[syn] <- "synonymous"
    
    # Missense or nonsense
    missOrNon <- which(!is.na(mvr[snp]$AAchange))
    for (i in 1:length(missOrNon)) {
      index <- missOrNon[i]
      aaChange <- unlist(strsplit(stri_replace(mvr[snp]$AAchange[index], " ", regex="\\d+"), " "))
      
      
      if ("*" %in% aaChange[2]) mvr[snp]$typeMut[index] <- "nonsense"
      else mvr[snp]$typeMut[index] <- "missense"
    }
    
  }
  
  ins <- grep("ins", names(mvr))
  if (length(ins) != 0) {
    
    insLength <- nchar(alt(mvr[ins])) - 1
    
    mvr[ins]$typeMut[which(insLength %% 3 == 0)] <- "insertion"
    mvr[ins]$typeMut[which(!(insLength %% 3 == 0))] <- "frameshift"
    
  }
  
  del <- grep("del", names(mvr))
  if (length(del) != 0) {
    
    delLength <- nchar(ref(mvr[del])) - 1
    
    mvr[del]$typeMut[which(delLength %% 3 == 0)] <- "deletion"
    mvr[del]$typeMut[which(!(delLength %% 3 == 0))] <- "frameshift"
  }
  
  return(mvr)
}
