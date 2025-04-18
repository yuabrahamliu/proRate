% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structure.R
\name{getexon}
\alias{getexon}
\title{Calculate exon coverages for specific gene regions}
\usage{
getexon(targetfile = NULL, genomename = "mm10", genenames = NULL)
}
\arguments{
\item{targetfile}{The directory of the file indicating the gene regions
whose exon coverages need to be calculated, not necessary to be exactly
the gene body region between TSS and TTS sites. Columns named as chr,
start, end, strand, and gene_id are required. If it is NULL, all the gene
body regions of UCSC known genes in the genome specified by the parameter
\code{genomename} will be analyzed.}

\item{genomename}{Specify the genome of the genes to be analyzed, when the
parameter \code{targetfile} is NULL.}

\item{genenames}{The symbols of genes whose exon coverages need to be
calculated. The regions indicated by the parameter \code{targetfile} or
\code{genomename} will be analyzed only if they belong to the genes
indicated by this parameter \code{genenames}.}
}
\value{
A data.frame with a column named exon indicating the exon coverages
of the gene regions.
}
\description{
Calculate exon coverages for specific genes regions.
}
