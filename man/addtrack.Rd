% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supps.R
\name{addtrack}
\alias{addtrack}
\title{Add gene coordinate track to a gene's rate inference result plot}
\usage{
addtrack(
  genedat,
  binplotdat,
  expandplotdat,
  genomename = "mm10",
  textsize = 13,
  titlesize = 15,
  face = "bold",
  method = NULL
)
}
\arguments{
\item{genedat}{The gene rate inference report data frame generated by the
function \code{calrate} or \code{mcalrate}. It can be extracted from the
"report" slot of these functions' result, which is a data frame, and what
needed here is the its sub-data-frame (with only one row) containing the
result information of the specific gene whose inference plots need to be
modified here.}

\item{binplotdat}{The bin-level rate inference plot data of a specific gene.
It can be extracted from the "binplots" slot of the results generated by
\code{calrate} or \code{mcalrate}. What needed here is the element data
contained in this slot using the gene's name as the element name.}

\item{expandplotdat}{The base-pair-level rate inference plot data of a gene.
It can be extracted from the "expandplots" slot of the results generated
by the function \code{calrate} or \code{mcalrate}. What needed here is the
element data contained in this slot using the gene's name as the element
name.}

\item{genomename}{Specify the genome of the specific gene whose inference
plots need to be modified here. It can be "mm10" for mouse or "hg38" for
human.}

\item{textsize}{The font size for the plot texts. Default is 13.}

\item{titlesize}{The font size for the plot titles. Default is 15.}

\item{face}{The font face for the plot texts. Default is "bold".}

\item{method}{The method used when inferring the gene's transcription rate
with the function \code{calrate} or \code{mcalrate}. Can be "LSS" for the
least sum of squares method, or "HMM" for the hidden Markov model method.}
}
\value{
The modified gene rate inference plots with the gene's coordinate
track added.
}
\description{
Add gene coordinate track to a gene's rate inference result plot generated
by the function \code{calrate} or \code{mcalrate}.
}
\examples{
library(proRate)

wt0file <- system.file("extdata", "wt0.bam", package = "proRate")
wt15file <- system.file("extdata", "wt15.bam", package = "proRate")

wtrates <- calrate(time1file = wt0file, 
                  time2file = wt15file, 
                  time = 15, 
                  strandmethod = 1, 
                  
                  genomename = "mm10", 
                  lencutoff = 40000, 
                  fpkmcutoff = 1, 
                  
                  threads = 1, 
                  
                  startshorten = 1000, 
                  endshorten = 1000, 
                  window_num = 40, 
                  
                  method = "LSS", 
                  pythonpath = NULL, 
                  
                  difftype = 1)

addtrack(genedat = subset(wtrates$report, gene_id == "Mamdc2"), 
        binplotdat = wtrates$binplots$Mamdc2, 
        expandplotdat = wtrates$expandplots$Mamdc2, 
        genomename = "mm10", 
        method = "LSS", 
        titlesize = 17, 
        textsize = 16, 
        face = "bold")



}
