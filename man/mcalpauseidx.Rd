% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pausing.R
\name{mcalpauseidx}
\alias{mcalpauseidx}
\title{Calculate Pause indeces for genes for multiple bam files}
\usage{
mcalpauseidx(
  bamfiles,
  genefile,
  genomename = NULL,
  tssradius = 1000,
  labels,
  strandmethod = 1,
  savegenenames = NULL,
  plotgenenames = TRUE,
  mergecases = FALSE,
  threads = 1,
  genelencutoff = NULL,
  fpkmcutoff = 1
)
}
\arguments{
\item{bamfiles}{The bam files need to calculate the Pause indeces for genes.
Should be a vector with elements as strings indicating the directores of
the files.}

\item{genefile}{The genes whose Pause indeces need to be calculated. If it
is NULL, all the genes in the genome specified by the parameter
\code{genomename} will be analyzed. If provied by the user, columns named
as chr, start, end, strand, and gene_id are required. All genes should be
longer than the maximum value of the parameters \code{tssradius} and
\code{genelencutoff}.}

\item{genomename}{Specify the genome of the genes to be analyzed, when the
parameter \code{genefile} is NULL.}

\item{tssradius}{A numeric value indicating the radius of the promoter
region centering around the TSS.}

\item{labels}{A vector with elements as strings to be included in the titles
of the specific gene plots to indicate the experimental conditions of the
bam files. There should be no replicated elements in this vector.}

\item{strandmethod}{The strand specific method used when preparing the
sequencing library, can be 1 for directional ligation method and 2 for
dUTP method. If the sample is sequenced using a single strand method, set
it as 0.}

\item{savegenenames}{For which genes their concrete FPM value for each bp
position need to be saved, or plotted.}

\item{plotgenenames}{Whether to plot the FPM value for each bp position for
the genes provided by the parameter \code{savegenenames}.}

\item{mergecases}{Whether merge the data of all bam files together to one,
and then use it to calculate gene Pause indeces and plot genes. Default is
FALSE.}

\item{threads}{Number of threads to do the parallelization. Default is 1.}

\item{genelencutoff}{The cutoff on gene length (bp). The default value is
NULL, but if it is set, only genes longer than this cutoff will be
considerred.}

\item{fpkmcutoff}{The cutoff on gene FPKM value. Only genes with an FPKM
value greater than the cutoff will be considerred. Default is 1.}
}
\value{
Will generate a list with several sub-lists recording the Pause
indeces of the genes, as well as plot the FPM value of each bp for the
specific genes indicated by the parameter \code{savegenenames}.
}
\description{
Calculate Pause indeces for genes for mulitple bam files (Pause index is the
ratio of transcription Polymerase II signal density near a gene promoter to
signal density within the gene body, such that higher Pause indices reflect
a greater enrichment of promoter-paused Polymerase II), and also plot the
FPM value of each bp for specific genes.
}
