% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/importVCF.R
\name{vcf2mvrl}
\alias{vcf2mvrl}
\title{Convert a VCF to an MVRangesList object}
\usage{
vcf2mvrl(vcf_filename, bam_filename = NULL, gzip = TRUE,
  biscuit = FALSE, verbose = TRUE)
}
\arguments{
\item{vcf_filename}{A filename/path to the corresponding VCF}

\item{bam_filename}{A filename/path to the corresponding BAM}

\item{gzip}{Is the input gzipped?}

\item{biscuit}{Is the input VCF from biscuit?}

\item{verbose}{logical, if set to TRUE, print all messages during progressing.}
}
\value{
an MVRangesList
}
\description{
Convert a VCF to an MVRangesList object
}
\examples{


}
