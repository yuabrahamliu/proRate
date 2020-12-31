#'Calculate GC content for specific gene regions
#'
#'Calculate GC content for specific genes regions.
#'
#'@param targetfile The directory of the file indicating the gene regions 
#'  whose GC contents need to be calculated, not necessary to be exactly the 
#'  gene body region between TSS and TTS sites. Columns named as chr, start, 
#'  end, strand, and gene_id are required. If it is NULL, all the gene body 
#'  regions of UCSC known genes in the genome specified by the parameter 
#'  \code{genomename} will be analyzed.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetfile} is NULL.
#'@param genenames The symbols of genes whose GC contents need to be 
#'  calculated. The default value is NULL. If it is set, for the regions 
#'  indicated by the parameter \code{targetfile} or \code{genomename}, they 
#'  will be analyzed only when they also belong to the genes indicated by this 
#'  parameter \code{genenames}. 
#'@return A data.frame with a column named GC indicating the GC contents of 
#'  the gene regions. 
#'@export
getgc <- function(targetfile, genomename = 'mm10', genenames = NULL){
  
  #library(GenomicAlignments)
  
  if(!is.null(targetfile)){
    
    if(file.exists(targetfile)){
      
      genes <- read.table(targetfile, sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE, quote = '',
                          check.names = FALSE)
      
      genes$start <- genes$start + 1
      
      if(!is.null(genenames)){
        genes <- subset(genes, gene_id %in% genenames)
      }
      
      geneframe <- genes
      
      genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      genes <- GenomicRanges::GRanges(seqnames = genes$seqnames,
                                      ranges = IRanges::IRanges(start = genes$start, end = genes$end),
                                      strand = genes$strand,
                                      gene_id = genes$gene_id)
      
      genes <- genegccont(origenes = genes, genomename = genomename)
      
      genes <- GenomicAlignments::sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- get(paste0(genomename, '.packagegenes'))
    genes <- GenomicRanges::GRanges(genes)
    
    if(!is.null(genenames)){
      genes <- subset(genes, gene_id %in% genenames)
    }
    
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- GenomeInfoDb::keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
                                         pruning.mode = 'coarse')
    
    genes$updis <- NULL
    genes$dndis <- NULL
    genes$exon <- NULL
    
  }
  
  geneframe <- as.data.frame(genes)
  
  for(i in 1:ncol(geneframe)){
    
    if(is.factor(geneframe[,i])){
      geneframe[,i] <- as.character(geneframe[,i])
    }
    
  }
  
  return(geneframe)
  
}



extractfeature <- function(genes = genes1,
                           feature = feature, radius = radius){
  
  starts <- genes@ranges@start
  ends <- starts + genes@ranges@width - 1
  strands <- as.character(genes@strand)
  
  if(feature == 'promoter'){
    centers <- starts
    centers[strands == '-'] <- ends[strands == '-']
  }else{
    centers <- ends
    centers[strands == '-'] <- starts[strands == '-']
  }
  
  featurestarts <- centers - radius
  featurestarts[featurestarts < 1] <- 1
  featureends <- centers + radius
  
  featureranges <- IRanges::IRanges(start = featurestarts, end = featureends)
  featuregrs <- GenomicRanges::GRanges(seqnames = genes@seqnames, ranges = featureranges,
                                       strand = strands, gene_id = genes$gene_id)
  
  return(featuregrs)
  
}

genekmer <- function(genes = genes, genomename = 'mm10', k = 6){
  
  #library(GenomicRanges)
  #library(Biostrings)
  
  geneends <- genes@ranges@start + genes@ranges@width - 1
  
  if(genomename == 'mm10'){
    #library("BSgenome.Mmusculus.UCSC.mm10")
    chrlimits <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus@seqinfo@seqlengths
    names(chrlimits) <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus@seqinfo@seqnames
    
    genelimits <- as.vector(chrlimits[as.character(genes@seqnames)])
    geneends[geneends > genelimits] <- genelimits[geneends > genelimits]
    genes <- GenomicRanges::GRanges(seqnames = genes@seqnames,
                                    ranges = IRanges::IRanges(start = genes@ranges@start,
                                                              end = geneends),
                                    strand = genes@strand,
                                    gene_id = genes$gene_id)
    
    seqs <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10::Mmusculus, genes)
  }else if(genomename == 'hg38'){
    #library("BSgenome.Hsapiens.UCSC.hg38")
    chrlimits <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens@seqinfo@seqlengths
    names(chrlimits) <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens@seqinfo@seqnames
    
    genelimits <- as.vector(chrlimits[as.character(genes@seqnames)])
    geneends[geneends > genelimits] <- genelimits[geneends > genelimits]
    genes <- GenomicRanges::GRanges(seqnames = genes@seqnames,
                                    ranges = IRanges::IRanges(start = genes@ranges@start,
                                                              end = geneends),
                                    strand = genes@strand,
                                    gene_id = genes$gene_id)
    
    seqs <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, genes)
  }
  
  counts_genes <- Biostrings::oligonucleotideFrequency(seqs, width = k, step = 1)
  counts_genes <- as.data.frame(counts_genes)
  freq_genes <- colSums(counts_genes)
  freq_genes <- freq_genes[order(-freq_genes)]
  
  return(freq_genes)
  
}

selectfisher <- function(freqvec1 = kmerres2, freqvec2 = kmerres1){
  
  A11 <- sum(freqvec1)
  A21 <- sum(freqvec2)
  
  fpvec <- c()
  kmers <- names(freqvec1)
  freqvec2 <- freqvec2[kmers]
  
  for(i in 1:length(kmers)){
    A12 <- freqvec1[kmers[i]]
    A22 <- freqvec2[kmers[i]]
    ftab <- matrix(c(A11, A12, A21, A22), ncol=2, nrow=2, byrow = TRUE)
    f <- fisher.test(ftab)
    fp <- f$p.value
    
    fpvec <- c(fpvec, fp)
    names(fpvec)[length(fpvec)] <- kmers[i]
  }
  
  fadjpvec <- p.adjust(fpvec, method = 'BH')
  
  kmerratio1 <- freqvec1/A11
  kmerratio2 <- freqvec2/A21
  
  kmerratio <- kmerratio1/kmerratio2
  
  res <- data.frame(kmer = kmers, ratio = kmerratio,
                    pval = fpvec, padj = fadjpvec,
                    stringsAsFactors = FALSE)
  res <- res[order(-res$ratio),]
  row.names(res) <- 1:nrow(res)
  
  return(res)
  
}



#'Compare the kmer difference between 2 sets of sequences
#'
#'Compare the kmer difference between 2 sets of sequences, reporting the ratio 
#'between the kmer frequencies of the 2 sequence sets, the p-value and 
#'adjusted p-value of the difference.
#'
#'@param targetfile The directory of the file indicating the gene regions 
#'  whose kmer need to be compared between the 2 sequence sets, not necessary 
#'  to be extactly the gene body region between TSS and TTS sites. Columns 
#'  named as chr, start, end, strand, and gene_id are required. If it is NULL, 
#'  all the gene regions defined together by the parameters \code{genomename}, 
#'  \code{feature} and \code{radius} will be analyzed.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetfile} is NULL.
#'@param k The length of the kmer to be analyzed. Default is 6, meaning 6-mer 
#'  will be analyzed. 
#'@param genes1 The symbols of genes in set 1 whose kmers need to be compared 
#'  with that of set 2. The regions indicated by the parameter 
#'  \code{targetfile} or \code{genomename} etc will be used for these genes 
#'  only if they belong to the genes indicated by this parameter 
#'  \code{genes1}. 
#'@param genes2 The symbols of genes in set 2 whose kmers need to be compared 
#'  with that of set 1. Similar to the parameter \code{genes1}.
#'@param feature If the parameter \code{targetfile} is NULL, while the 
#'  parameter \code{genomename} is defined. This parameter \code{feature} can 
#'  be used to further select regions from the genes indicated by 
#'  \code{genomename}. Can choose from 'promoter', 'end', and 'genebody'. If 
#'  it is 'promoter' or 'end', another parameter \code{radius} is needed to 
#'  define the radius of the promoter or end region centering around the TSS 
#'  or TTS site.
#'@param radius A numberic value needed to define the radius length of the 
#'  gene promoter or end region if the parameter \code{feature} is set as 
#'  'promoter' or 'end'.
#'@return A data.frame indicating the kmer frquency ratios between sequence 
#'  set 1 and set 2, as well as p-values (calculated with Fisher's test) and 
#'  adjusted p-values (adjusted using the Benjamini & Hochberg method).
#'@export
getkmer <- function(targetfile, genomename = 'mm10', k = 6,
                    genes1,
                    genes2,
                    feature,
                    radius = 1000){
  
  #library(GenomicAlignments)
  
  if(!is.null(targetfile)){
    
    if(file.exists(targetfile)){
      
      genes <- read.table(targetfile, sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE, quote = '',
                          check.names = FALSE)
      
      genes$start <- genes$start + 1
      
      geneframe <- genes
      
      genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      genes <- GenomicRanges::GRanges(seqnames = genes$seqnames,
                                      ranges = IRanges::IRanges(start = genes$start, end = genes$end),
                                      strand = genes$strand,
                                      gene_id = genes$gene_id)
      
      genes <- GenomicAlignments::sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- get(paste0(genomename, '.packagegenes'))
    genes <- GenomicRanges::GRanges(genes)
    
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- GenomeInfoDb::keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
                                         pruning.mode = 'coarse')
    
    genes$updis <- NULL
    genes$dndis <- NULL
    genes$GC <- NULL
    genes$exon <- NULL
    
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    if(feature != 'genebody'){
      genes <- extractfeature(genes = genes, feature = feature, radius = radius)
    }else{
      genes <- genes
    }
    
  }
  
  genes1 <- subset(genes, gene_id %in% genes1)
  genes2 <- subset(genes, gene_id %in% genes2)
  
  genes1 <- GenomicAlignments::sort(genes1)
  genes2 <- GenomicAlignments::sort(genes2)
  names(genes1) <- 1:length(genes1)
  names(genes2) <- 1:length(genes2)
  
  kmerres1 <- genekmer(genes = genes1, genomename = genomename, k = k)
  kmerres2 <- genekmer(genes = genes2, genomename = genomename, k = k)
  
  kmerres <- selectfisher(freqvec1 = kmerres2, freqvec2 = kmerres1)
  
  return(kmerres)
  
}



#'Calculate exon coverages for specific gene regions
#'
#'Calculate exon coverages for specific genes regions.
#'
#'@param targetfile The directory of the file indicating the gene regions 
#'  whose exon coverages need to be calculated, not necessary to be exactly 
#'  the gene body region between TSS and TTS sites. Columns named as chr, 
#'  start, end, strand, and gene_id are required. If it is NULL, all the gene 
#'  body regions of UCSC known genes in the genome specified by the parameter 
#'  \code{genomename} will be analyzed.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetfile} is NULL.
#'@param genenames The symbols of genes whose exon coverages need to be 
#'  calculated. The regions indicated by the parameter \code{targetfile} or 
#'  \code{genomename} will be analyzed only if they belong to the genes 
#'  indicated by this parameter \code{genenames}. 
#'@return A data.frame with a column named exon indicating the exon coverages 
#'  of the gene regions. 
#'@export
getexon <- function(targetfile, genomename = 'mm10', genenames = NULL){
  
  #library(GenomicAlignments)
  
  if(!is.null(targetfile)){
    
    if(file.exists(targetfile)){
      
      genes <- read.table(targetfile, sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE, quote = '',
                          check.names = FALSE)
      
      genes$start <- genes$start + 1
      
      if(!is.null(genenames)){
        genes <- subset(genes, gene_id %in% genenames)
      }
      
      geneframe <- genes
      
      genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      genes <- GenomicRanges::GRanges(seqnames = genes$seqnames,
                                      ranges = IRanges::IRanges(start = genes$start, end = genes$end),
                                      strand = genes$strand,
                                      gene_id = genes$gene_id)
      
      genes <- geneexonfreq(origenes = genes, genomename = genomename)
      
      genes <- GenomicAlignments::sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- get(paste0(genomename, '.packagegenes'))
    genes <- GenomicRanges::GRanges(genes)
    
    if(!is.null(genenames)){
      genes <- subset(genes, gene_id %in% genenames)
    }
    
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- GenomeInfoDb::keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
                                         pruning.mode = 'coarse')
    
    genes$updis <- NULL
    genes$dndis <- NULL
    genes$GC <- NULL
    
  }
  
  geneframe <- as.data.frame(genes)
  
  for(i in 1:ncol(geneframe)){
    
    if(is.factor(geneframe[,i])){
      geneframe[,i] <- as.character(geneframe[,i])
    }
    
  }
  
  return(geneframe)
  
}


