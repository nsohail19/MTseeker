#' Plot an empty circos plot of the mitochondrial genome
#'
#' @param mvr      MVRanges
#'
#' @return        Image of mitochondrial genome
#'
#'
#' @examples 
#' 
#' library(MTseekerData)
#' BAMdir <- system.file("extdata", "BAMs", package="MTseekerData")
#' BAMs <- list.files(BAMdir, pattern="bam$")
#' BAMs <- paste0(BAMdir, "/", BAMs)
#' 
#' pu <- pileupMT(BAMs[1], ref="rCRS)
#' genMTcircos(pu)
#' 
#' 
#'
#' @export

genMTcircos <- function(mvr) {
  
  circos.clear()
  
  anno <- initMTcircos(mvr)
  genesMTcircos(mvr, anno, legends=T)
  
}

# Initialize the circos plots
initMTcircos <- function(x) {

  # Human
  if (unique(genome(x)) == "rCRS") {
    data(mtAnno.rCRS)
    anno <- mtAnno #.rCRS
    refWidth <- 16569
  }
  
  # Mouse
  else if (unique(genome(x)) == "NC_005089") {
    anno <- readRDS("~/Documents/pileupTesting/NC_005089genome/MTmouseAnno.rds")
    refWidth <- 16299
  }
  
  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  circos.par("clock.wise"=FALSE, "start.degree"=90, "gap.degree"=0, 
             "track.margin"=c(0.005, 0.005), "cell.padding"=c(0.005,0,0.005,0), 
             "points.overflow.warning"=FALSE)
  circos.genomicInitialize(data=dat, plotType=NULL, major.by=refWidth)
  
  return(anno)
}

# Actually plot the rings with the gene names
genesMTcircos <- function(x, anno, legends=F) {

  dat <- data.frame(name=names(anno), start=start(anno), end=end(anno))
  
  # Needs to go inside the function for some reason
  pfun <- function(x, y) {
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    gr <- anno[CELL_META$sector.index]
    ytop <- .height(gr) * ifelse(strand(gr) == "+", 1, 0)
    ybot <- .height(gr) * ifelse(strand(gr) == "-", -1, 0)
    lab <- ifelse(CELL_META$sector.index == "DLP", "CR", CELL_META$sector.index)
    circos.rect(xlim[1], ybot, xlim[2], ytop, col=gr$itemRgb)
    if (gr$region %in% c("rRNA", "coding", "D-loop") & gr$name != "HVR3") {
      circos.text(mean(xlim), .textloc(gr), lab, col="black", cex=.textcex(gr),
                  font=.textbold(gr), facing="clockwise", niceFacing=TRUE)
    }
  }
  
  circos.track(panel.fun=pfun, ylim=c(-1,1), track.height=0.5, 
               track.margin=c(0,0), bg.border=NA)
  
  if (legends && genome(x) == "rCRS") {
    
    # Get the colors to correspond which region they belong to
    # Only support rCRS 
    colDF <- as.data.frame(unique(anno$itemRgb), stringsAsFactors=F)
    names(colDF) <- "col"
    
    colDF$label <- NA_character_
    colDF$label[1] <- "Control Region"
    colDF$label[2] <- "D Loop"
    colDF$label[3] <- "tRNA"
    colDF$label[4] <- "rRNA"
    colDF$label[5] <- "Complex I"
    colDF$label[6] <- "Complex IV"
    colDF$label[7] <- "Complex V"
    colDF$label[8] <- "Complex III"
    
    legend("topright", title="Regions",
           legend=colDF$label, col=colDF$col, pch=15, cex=0.8)
  }
  
  res <- list(anno=dat, pfun=pfun)
  return(res)
}

# helper fn
.height <- function(gr) ifelse(gr$region == "tRNA", 0.5, 1)

# helper fn
.halfheight <- function(gr) gr$region == "tRNA"

# helper fn
.textloc <- function(gr) {
  ifelse(.halfheight(gr), .25, .5) * ifelse(strand(gr) == "+", 1, -1)
}

# helper fn
.textbold <- function(gr) ifelse(gr$region %in% c("coding", "rRNA"), 3, 1)

# helper fn
.textcex <- function(gr) { 
  ifelse(gr$name %in% c("HVR1","HVR2","MT-ATP8"), .65,
         ifelse(gr$name %in% c("MT-ND3","MT-ND4L","MT-ND6","MT-CO2","MT-CO3",
                               "MT-RNR1","MT-RNR2","MT-ATP8","MT-ATP6"),.8,.95))
}

# helper fn
.colorCode <- function(x, darken=TRUE, howMuch=1.25) { 
  data("mtAnno.rCRS", package="MTseeker")
  color <- mtAnno[x]$itemRgb
  if (darken) color <- .darken(color, howMuch=howMuch)
  return(color)
}

# helper fn
.darken <- function(hex, howMuch=1.25) {
  rgb(t(col2rgb(hex)/howMuch), maxColorValue=255)
}

# helper fn
.newsprint <- colorRamp2(c(0, 1), c("#FFFFFF", "#000000"))

# helper fn
.bloody <- colorRamp2(c(0, 1), c("#FFFFFF", "#880000"))

# helper fn
.blurple <- colorRamp2(c(0, 1), c("#FFFFFF", "#FF00FF"))

# helper fn
.viridis <- colorRamp2(seq(0, 1, by = 0.1), viridis(11))

# helper fn
.plasma <- colorRamp2(seq(0, 1, by = 0.1), plasma(11))

# helper fn; should probably use viridis instead 
.jet <- colorRamp2(seq(0, 1, 0.125),
                   c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
