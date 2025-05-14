#rm(list = ls())

#gc()

#wkdir <- 'C:/Users/yuabr/Desktop/Transfer/codetransfer/proRate/proRate_V2/'

#setwd(wkdir)

#Import data####

#library(proRate)

#wt0file <- system.file('extdata', 'wt0.bam', package = 'proRate')
#wt15file <- system.file('extdata', 'wt15.bam', package = 'proRate')
#ko0file <- system.file('extdata', 'ko0.bam', package = 'proRate')
#ko15file <- system.file('extdata', 'ko15.bam', package = 'proRate')

#targetfile <- system.file('extdata', 'targetgenes.txt', package = 'proRate')

#fp0file <- './testdatasets/fp_00_bamAligned.sortedByCoord.out.qc.chr19.bam'
#fp12file <- './testdatasets/fp_12_bamAligned.sortedByCoord.out.qc.chr19.bam'


#Functions####

#Some parallel computing going on in the background that is not getting cleaned up
#fully between runs can cause Error in summary.connection(connection) : invalid
#connection. The following function is needed to be called to fix this error.
unregister_dopar <- function(){
  
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  
}

removeoverlap <- function(coords = genecoords, ignorestrand = TRUE){
  
  dis <- GenomicRanges::distanceToNearest(coords, ignore.strand = ignorestrand)
  disvec <- dis@elementMetadata$distance
  
  if(sum(disvec == 0)  == 0){
    return(coords)
  }else{
    newgenecoords <- coords[disvec != 0]
    removeoverlap(coords = newgenecoords, ignorestrand = ignorestrand)
  }
  
}

geneintervalscreen <- function(origenes = genes, ignorestrand = TRUE){
  
  #library(GenomicRanges)
  
  if(ignorestrand == TRUE){
    ups <- GenomicRanges::follow(origenes, BiocGenerics::unstrand(origenes))
    dns <- GenomicRanges::precede(origenes, BiocGenerics::unstrand(origenes))
  }else{
    ups <- GenomicRanges::follow(origenes, origenes)
    dns <- GenomicRanges::precede(origenes, origenes)
  }
  #follow returns the nearest upstream region of the query and considers the strandness of the query.
  #If the query is '+' stranded, its upstream region returned will before its start (also on '+' strand,
  #or on either '+' or '-' strand, if the strand information of the subject is removed with unstrand())
  #If the query is '-' stranded, its upstream region returned will after its end (also on '-' strand,
  #or on either '-' or '+' strand, if the strand information of the subject is removed with unstrand())
  #precede returns the nearest downstream region of the query and considers the strandness of the query.
  #If the query is '+' stranded, its downstream region returned will after its end
  #If the query is '-' stranded, its downstream region returned will before its start
  
  strandups <- ups
  stranddns <- dns
  
  strandupdis <- GenomicFeatures::distance(x = origenes[!is.na(strandups)], y = origenes[strandups[!is.na(strandups)]],
                                           ignore.strand = ignorestrand)
  stranddndis <- GenomicFeatures::distance(x = origenes[!is.na(stranddns)], y = origenes[stranddns[!is.na(stranddns)]],
                                           ignore.strand = ignorestrand)
  #distance returns the distance between the closest 2 points of the upstream and the downstream regions
  #(the 2 points are not included in the distance), no matter both of them are on '+' strand or '-' strand.
  #However, if x and y are on different strands, it will return NA, unless to set ignore.strand = TRUE,
  #and it will treat them as on the same strand
  
  
  
  strandupdises <- rep(NA, length(origenes))
  stranddndises <- rep(NA, length(origenes))
  strandupdises[!is.na(strandups)] <- strandupdis
  stranddndises[!is.na(stranddns)] <- stranddndis
  origenes$updis <- strandupdises
  origenes$dndis <- stranddndises
  
  
  
  otherstrandupdis <- origenes[is.na(strandups)]
  otherstranddndis <- origenes[is.na(stranddns)]
  
  outlierdis <- function(targets = otherstrandupdis, dir = 'up'){
    
    chrlens <- targets@seqinfo@seqlengths
    names(chrlens) <- targets@seqinfo@seqnames
    
    targetchrlens <- chrlens[as.character(targets@seqnames)]
    disstart <- targets@ranges@start
    disend <- targetchrlens - (targets@ranges@start + targets@ranges@width - 1)
    disend <- as.vector(disend)
    
    dises <- matrix(c(disstart, disend), nrow = 2, byrow = TRUE)
    dises <- as.vector(dises)
    
    if(dir == 'up'){
      dises <- dises[rep(c('+', '-'), length(targets)) == rep(as.character(targets@strand), each = 2)]
      targets$updis <- dises
    }else{
      dises <- dises[rep(c('-', '+'), length(targets)) == rep(as.character(targets@strand), each = 2)]
      targets$dndis <- dises
    }
    
    return(targets)
    
  }
  
  otherstrandupdis <- outlierdis(targets = otherstrandupdis, dir = 'up')
  otherstranddndis <- outlierdis(targets = otherstranddndis, dir = 'dn')
  
  
  origenes$updis[is.na(origenes$updis)] <- otherstrandupdis$updis
  origenes$dndis[is.na(origenes$dndis)] <- otherstranddndis$dndis
  
  return(origenes)
  
}

genegccont <- function(origenes = genes, genomename = 'mm10'){
  
  #library(GenomicRanges)
  #library(Biostrings)
  
  if(genomename == 'mm10'){
    #library("BSgenome.Mmusculus.UCSC.mm10")
    seqs <- Biostrings::getSeq(BSgenome.Mmusculus.UCSC.mm10::Mmusculus, origenes)
  }else if(genomename == 'hg38'){
    #library("BSgenome.Hsapiens.UCSC.hg38")
    seqs <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, origenes)
  }
  
  #counts_genes <- oligonucleotideFrequency(seqs, width = 6, step = 1)
  #counts_genes <- as.data.frame(counts_genes)
  #freq_genes <- colSums(counts_genes)
  #freq_genes <- freq_genes[order(-freq_genes)]
  
  seqlist <- unlist(strsplit(toString(seqs), split = ', '))
  seqframe <- data.frame(genesym = origenes$gene_id,
                         seq = seqlist, stringsAsFactors = FALSE)
  
  #library(seqinr)
  calGC <- function(line){
    seqstring <- as.vector(unlist(line[2]))
    seqvector <- strsplit(seqstring, split = '')[[1]]
    gccontent <- seqinr::GC(seqvector)
    return(gccontent)
  }
  
  genegc <- apply(seqframe, 1, calGC)
  genegc <- data.frame(genesym = seqframe$genesym, seqgc = genegc, stringsAsFactors = FALSE)
  
  origenes$GC <- genegc$seqgc
  
  return(origenes)
  
}

calexonfreq <- function(genelineidx, origenes, exoncoords){
  
  geneline <- origenes[genelineidx]
  geneexons <- GenomicRanges::intersect(exoncoords, geneline)
  
  if(length(geneexons) == 0){
    return(0)
  }else{
    
    geneexons <- GenomicRanges::reduce(geneexons)
    geneexons$gene_id <- geneline$gene_id
    exonlens <- sum(geneexons@ranges@width)
    genelen <- geneline@ranges@width
    exonfreq <- exonlens/genelen
    
    return(exonfreq)
    
  }
  
}

geneexonfreq <- function(origenes = genes, genomename = 'mm10'){
  
  #origenes <- get(paste0(genomename, '.packagegenes'))[1:100]
  
  
  #library(AnnotationDbi)
  
  if(genomename == 'hg38'){
    
    #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    #library(org.Hs.eg.db)
    
    exoncoords <- GenomicFeatures::exons(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
    
  }else if(genomename == 'mm10'){
    
    #library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    #library(org.Mm.eg.db)
    
    exoncoords <- GenomicFeatures::exons(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene)
    
  }
  
  geneexonfreqs <- sapply(X = 1:length(origenes), FUN = calexonfreq, 
                          origenes = origenes, exoncoords = exoncoords)
  
  origenes$exon <- geneexonfreqs
  
  return(origenes)
  
}

getfpkm <- function(features = genes, 
                    reads = reads1,
                    readsdepth = reads1depth, 
                    singleend = FALSE){
  
  
  fcounts <- GenomicAlignments::summarizeOverlaps(features = features, reads = reads,
                                                  mode = 'IntersectionNotEmpty',
                                                  singleEnd = singleend, 
                                                  ignore.strand = FALSE)
  #Add singleEnd = FALSE for pair-end
  #Add ignore.strand = FALSE for strand specific
  
  fcounts <- as.data.frame(SummarizedExperiment::assays(fcounts)$counts)
  #assays(x, ...)
  #The SummarizedExperiment class is a matrix-like container where rows represent 
  #features of interest (e.g. genes, transcripts, exons, etc...) and columns represent 
  #samples (with sample data summarized as a DataFrame). A SummarizedExperiment object 
  #contains one or more assays, each represented by a matrix-like object of numeric
  #or other mode.
  #x A SummarizedExperiment object
  
  featurewidth <- features@ranges@width
  
  factorM <- readsdepth/10^6
  factorK <- featurewidth/10^3
  
  fpkms <- fcounts/factorK/factorM
  features$fpkm <- fpkms$reads
  
  return(features)
  
}

#'Generate a GRanges object of UCSC known genes in specific genome
#'
#'Generate a GRanges object of UCSC known genes in specific genome.
#'
#'@param genomename The name of the genome, can be "mm10" for mouse or "hg38" 
#'  for human.
#'@param save Whether save the result GRanges as a file in work directory. 
#'  Default is TRUE.
#'@param calgenegccont Whether need to calculate the GC content for each gene 
#'  and include it into the result GRanges object. Default is TRUE.
#'@param calgeneexonfreq Whether calculate the exon frequency for each gene 
#'  and include it into the result GRanges object. Default is TRUE.
#'@return A GRanges object of UCSC known genes, with several gene information 
#'  provided, including seqnames, ranges, strands, gene symbols, distance to 
#'  the nearest upstream gene, distance to the nearest downstream gene, etc. 
#'  Genes overlapping with each other will not be included.
#'@export
makegenelist <- function(genomename = 'mm10', 
                         save = TRUE, 
                         
                         calgenegccont = TRUE, 
                         calgeneexonfreq = TRUE){
  
  #library(AnnotationDbi)
  
  if(genomename == 'hg38'){
    
    #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    #library(org.Hs.eg.db)
    
    genecoords <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
    #Select genes with unique gene symbols
    gene_ids <- genecoords$gene_id
    genesyms <- AnnotationDbi::select(x = org.Hs.eg.db::org.Hs.eg.db, keys = gene_ids,
                                      columns = 'SYMBOL', keytype = 'ENTREZID')
    
    targetchrs <- paste0('chr', c(seq(1, 22, 1), 'X', 'Y', 'M'))
    
  }else if(genomename == 'mm10'){
    
    #library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    #library(org.Mm.eg.db)
    
    genecoords <- GenomicFeatures::genes(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene)
    #Select genes with unique gene symbols
    gene_ids <- genecoords$gene_id
    genesyms <- AnnotationDbi::select(x = org.Mm.eg.db::org.Mm.eg.db, keys = gene_ids,
                                      columns = 'SYMBOL', keytype = 'ENTREZID')
    
    targetchrs <- paste0('chr', c(seq(1, 19, 1), 'X', 'Y', 'M'))
    
  }
  
  
  genesyms <- genesyms[complete.cases(genesyms), , drop = FALSE]
  names(genesyms) <- c('entrez', 'gene')
  dupgenes <- unique(genesyms$gene[duplicated(genesyms$gene)])
  genesyms <- subset(genesyms, !(gene %in% dupgenes))
  
  genecoords <- genecoords[gene_ids %in% genesyms$entrez]
  genecoords$gene_id <- genesyms$gene
  
  
  #Select genes only on chr1 to chrM
  genecoords <- GenomeInfoDb::keepSeqlevels(x = genecoords, value = targetchrs, pruning.mode = 'coarse')
  
  #Remove overlapped genes, no matter the strands are different or not###
  genecoords <- removeoverlap(coords = genecoords, ignorestrand = TRUE)
  
  #Add the information of the upstream distance and downstream distance to the nearest genes,
  #no matter the strands are different or not###
  genecoords <- geneintervalscreen(origenes = genecoords, ignorestrand = TRUE)
  
  
  genecoords <- GenomeInfoDb::keepSeqlevels(x = genecoords, value = unique(as.character(genecoords@seqnames)),
                                            pruning.mode = 'coarse')
  
  if(calgenegccont == TRUE){
    genecoords <- genegccont(origenes = genecoords, genomename = genomename)
  }
  
  if(calgeneexonfreq == TRUE){
    genecoords <- geneexonfreq(origenes = genecoords, genomename = genomename)
  }
  
  
  genecoords <- GenomicAlignments::sort(genecoords)
  genes <- genecoords
  names(genes) <- 1:length(genes)
  
  if(save == TRUE){
    saveRDS(genes, paste0(genomename, '.packagegenes.rds'))
  }
  
  return(genes)
  
}

#mm10genes <- makegenelist(genomename = 'mm10')

#hg38genes <- makegenelist(genomename = 'hg38')



getlastexon <- function(genelineidx = 1, origenes, exoncoords){
  
  geneline <- origenes[genelineidx]
  geneexons <- GenomicRanges::intersect(exoncoords, geneline)
  
  if(length(geneexons) == 0){
    
    lastexon <- NULL
    
  }else{
    
    if(as.vector(geneline@strand) == '+'){
      
      lastexon <- geneexons[length(geneexons)]
      lastexon$gene_id <- geneline$gene_id
      
    }else if(as.vector(geneline@strand) == '-'){
      
      lastexon <- geneexons[1]
      lastexon$gene_id <- geneline$gene_id
      
    }else{
      
      lastexon <- NULL
      
    }
    
  }
  
  return(lastexon)
  
}

makeutrlist <- function(genes, 
                        genomename = 'mm10'){
  
  #library(AnnotationDbi)
  
  if(genomename == 'hg38'){
    
    #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    #library(org.Hs.eg.db)
    
    exoncoords <- GenomicFeatures::exons(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
    
  }else if(genomename == 'mm10'){
    
    #library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    #library(org.Mm.eg.db)
    
    exoncoords <- GenomicFeatures::exons(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene)
    
  }
  
  lastexons <- sapply(X = 1:length(genes), FUN = getlastexon, origenes = genes, exoncoords = exoncoords)
  
  lastexons <- do.call(c, lastexons)
  
  return(lastexons)
  
}



parsereadsfile <- function(readsfile = time1file, strandmethod = 0){
  
  #library(GenomicAlignments)
  
  if(is.character(readsfile)){
    if(file.exists(readsfile)){
      
      if(strandmethod %in% c(0, 1, 2)){
        
        flags <- Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, 
                                        isDuplicate = FALSE, 
                                        isSecondaryAlignment = FALSE, 
                                        
                                        isProperPair = TRUE)
        
        params <- Rsamtools::ScanBamParam(flag = flags)
        
        reads_ori <- GenomicAlignments::readGAlignmentPairs(readsfile, 
                                                            strandMode = strandmethod, 
                                                            
                                                            param = params)
        
      }else{
        
        flags <- Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, 
                                        isDuplicate = FALSE, 
                                        isSecondaryAlignment = FALSE)
        
        params <- Rsamtools::ScanBamParam(flag = flags)
        
        reads_ori <- GenomicAlignments::readGAlignments(readsfile, 
                                                        
                                                        param = params)
      }
      reads <- as(reads_ori, 'GRanges')
      
    }else{
      
      reads <- GenomicRanges::GRanges()
      
      return(reads)
      
    }
    
  }else if(class(readsfile)[1] == 'GAlignmentPairs'){
    
    reads <- as(readsfile, 'GRanges')
    
  }else if(class(readsfile)[1] == 'GAlignments'){
    
    reads <- as(readsfile, 'GRanges')
    
  }else if(class(readsfile)[1] == 'GRanges'){
    
    reads <- readsfile
    
  }
  
  if(!('UCSC' %in% GenomeInfoDb::seqlevelsStyle(reads))){
    
    GenomeInfoDb::seqlevelsStyle(reads) <- 'UCSC'
    
  }
  
  if(length(reads) == 0){
    
    return(reads)
    
  }
  
  reads <- reads[grepl(pattern = 'chr[0-9].*|chrX|chrY', 
                       x = GenomeInfoDb::seqnames(reads))]
  
  #readscutoff <- max(301, min(quantile(reads@ranges@width, 0.95), 1000))
  
  #reads <- reads[reads@ranges@width < readscutoff]
  
  return(reads)
  
}



getreadslist <- function(timefiles = time1files, 
                         mergefiles = mergerefs, 
                         strandmode = strandmethod){
  
  readses <- list()
  
  i <- 1
  for(i in 1:length(timefiles)){
    timefile <- timefiles[i]
    
    if(mergefiles == TRUE){
      
      reads <- parsereadsfile(timefile, strandmethod = strandmode)
      
    }else{
      
      reads <- timefile
      
    }
    
    readses[[i]] <- reads
    
  }
  
  if(mergefiles == TRUE){
    readses <- do.call(c, readses)
    readses <- list(readses)
    
  }
  
  
  return(readses)
  
  
  
}



dirscan <- function(i = 8699, 
                    ext, 
                    
                    reads1, 
                    reads2, 
                    reads1depth, 
                    reads2depth){
  
  ext <- ext[i]
  
  chr <- GenomicAlignments::seqlevelsInUse(ext)
  GenomeInfoDb::seqlevels(ext) <- chr
  
  strandsign <- as.character(ext@strand@values)
  
  if(strandsign == '+'){
    
    start <- ext@ranges@start
    end <- start + ext@ranges@width - 1
    
    extreads1 <- IRanges::subsetByOverlaps(x = reads1, 
                                           ranges = ext)
    
    extreads2 <- IRanges::subsetByOverlaps(x = reads2, 
                                           ranges = ext)
    
    GenomeInfoDb::seqlevels(extreads1) <- chr
    GenomeInfoDb::seqlevels(extreads2) <- chr
    
    extreadscounts1 <- GenomicAlignments::coverage(extreads1, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    
    extreadscounts2 <- GenomicAlignments::coverage(extreads2, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    
    extreadscounts1 <- extreadscounts1[[1]]
    extreadscounts2 <- extreadscounts2[[1]]
    
  }else{
    
    start <- ext@ranges@start
    end <- start + ext@ranges@width - 1
    
    extreads1 <- IRanges::subsetByOverlaps(x = reads1, 
                                           ranges = ext)
    
    extreads2 <- IRanges::subsetByOverlaps(x = reads2, 
                                           ranges = ext)
    
    GenomeInfoDb::seqlevels(extreads1) <- chr
    GenomeInfoDb::seqlevels(extreads2) <- chr
    
    extreadscounts1 <- GenomicAlignments::coverage(extreads1, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    #Calculate the coverage per base, like wig file, per chrom
    
    extreadscounts2 <- GenomicAlignments::coverage(extreads2, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    
    extreadscounts1[[1]] <- rev(extreadscounts1[[1]])
    extreadscounts2[[1]] <- rev(extreadscounts2[[1]])
    
    extreadscounts1 <- extreadscounts1[[1]]
    extreadscounts2 <- extreadscounts2[[1]]
    
  }
  
  extdepth1 <- sum(extreadscounts1)
  extdepth2 <- sum(extreadscounts2)
  
  if(extdepth1 == 0 | extdepth2 == 0){
    return(NULL)
  }
  
  extreadscounts1_adj <- extreadscounts1/(reads1depth/10^6)
  extreadscounts2_adj <- extreadscounts2/(reads2depth/10^6)
  
  extdepth1_adj <- sum(extreadscounts1_adj)
  extdepth2_adj <- sum(extreadscounts2_adj)
  
  weight1 <- extdepth1_adj/mean(c(extdepth1_adj, extdepth2_adj))
  weight2 <- extdepth2_adj/mean(c(extdepth1_adj, extdepth2_adj))
  
  extreadscounts1_adj <- extreadscounts1_adj/weight1
  extreadscounts2_adj <- extreadscounts2_adj/weight2
  
  #Add pseudo count!!!
  #extreadscounts1_adj <- extreadscounts1_adj + 10^-5/weight1
  #extreadscounts2_adj <- extreadscounts2_adj + 10^-5/weight2
  
  #Add pseudo count 1!!!
  extreadscounts1_adj <- extreadscounts1_adj + 1/weight1
  extreadscounts2_adj <- extreadscounts2_adj + 1/weight2
  
  extreadscounts_adj <- extreadscounts1_adj + extreadscounts2_adj
  
  scanextreads <- as.vector(extreadscounts_adj)
  
  windowsizes <- seq(from = 50, 
                     to = ceiling(length(scanextreads)/2/50)*50, 
                     by = 50)
  
  windowdiffs <- c()
  
  windowmeanses <- list()
  
  for(m in 1:length(windowsizes)){
    
    windowsize <- windowsizes[m]
    
    windowmeans <- c()
    
    for(n in 1:length(scanextreads)){
      
      windowmean <- mean(scanextreads[n:min((windowsize + n - 1), 
                                            length(scanextreads))])
      windowmeans <- c(windowmeans, windowmean)
      
      cat(paste0('m = ', m, '; n = ', n, '\n'))
      
    }
    
    windowmeanses[[m]]  <- windowmeans
    
    windowdiff <- windowmeans[1] - min(windowmeans)
    
    windowdiffs <- c(windowdiffs, windowdiff)
    
  }
  
  windowidx <- which.max(windowdiffs)
  
  endpoint <- which.min(windowmeanses[[windowidx]]) - 1
  
  if(endpoint == 0){
    
    return(NULL)
    
  }
  
  if(strandsign == '+'){
    
    endpoint <- start + endpoint - 1
    
    endext <- GenomicRanges::GRanges(seqnames = chr, 
                                     ranges = IRanges::IRanges(start = start, 
                                                               end = endpoint), 
                                     strand = strandsign)
    
    endext@elementMetadata <- ext@elementMetadata
    
  }else{
    
    endpoint <- end - endpoint + 1
    
    endext <- GenomicRanges::GRanges(seqnames = chr, 
                                     ranges = IRanges::IRanges(start = endpoint, 
                                                               end = end), 
                                     strand = strandsign)
    
    endext@elementMetadata <- ext@elementMetadata
    
  }
  
  
  windowsize <- windowsizes[windowidx]
  
  res <- endext
  
  return(res)
  
}

mixscan <- function(i = 7302, 
                    mixext, 
                    extpair, 
                    
                    reads1, 
                    reads2, 
                    reads1depth, 
                    reads2depth){
  
  ext1 <- mixext[i]
  extpair <- subset(extpair, mixexts1 == ext1$gene_id)
  ext2 <- mixext[mixext$gene_id == extpair$mixexts2]
  
  chr <- GenomicAlignments::seqlevelsInUse(ext1)
  GenomeInfoDb::seqlevels(ext1) <- chr
  GenomeInfoDb::seqlevels(ext2) <- chr
  
  strandsign1 <- as.character(ext1@strand@values)
  strandsign2 <- as.character(ext2@strand@values)
  
  strandsigns <- c(strandsign1, strandsign2)
  exts <- c(ext1, ext2)
  names(exts) <- c(1, 2)
  
  mergedext <- GenomicRanges::reduce(exts, ignore.strand = TRUE)
  
  if(strandsign1 == '+'){
    
    start <- mergedext@ranges@start
    end <- start + mergedext@ranges@width - 1
    
    extreads1 <- IRanges::subsetByOverlaps(x = reads1, 
                                           ranges = mergedext)
    
    extreads2 <- IRanges::subsetByOverlaps(x = reads2, 
                                           ranges = mergedext)
    
    GenomeInfoDb::seqlevels(extreads1) <- chr
    GenomeInfoDb::seqlevels(extreads2) <- chr
    
    extreadscounts1 <- GenomicAlignments::coverage(extreads1, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    
    extreadscounts2 <- GenomicAlignments::coverage(extreads2, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    
    extreadscounts1 <- extreadscounts1[[1]]
    extreadscounts2 <- extreadscounts2[[1]]
    
  }else{
    
    start <- mergedext@ranges@start
    end <- start + mergedext@ranges@width - 1
    
    extreads1 <- IRanges::subsetByOverlaps(x = reads1, 
                                           ranges = mergedext)
    
    extreads2 <- IRanges::subsetByOverlaps(x = reads2, 
                                           ranges = mergedext)
    
    GenomeInfoDb::seqlevels(extreads1) <- chr
    GenomeInfoDb::seqlevels(extreads2) <- chr
    
    extreadscounts1 <- GenomicAlignments::coverage(extreads1, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    #Calculate the coverage per base, like wig file, per chrom
    
    extreadscounts2 <- GenomicAlignments::coverage(extreads2, 
                                                   shift = -(start - 1), 
                                                   width = (end - start + 1))
    
    extreadscounts1[[1]] <- rev(extreadscounts1[[1]])
    extreadscounts2[[1]] <- rev(extreadscounts2[[1]])
    
    extreadscounts1 <- extreadscounts1[[1]]
    extreadscounts2 <- extreadscounts2[[1]]
    
  }
  
  extdepth1 <- sum(extreadscounts1)
  extdepth2 <- sum(extreadscounts2)
  
  if(extdepth1 == 0 | extdepth2 == 0){
    
    newextgenes <- exts
    
    scanextreadses <- list()
    
    res <- list()
    res$extgenes <- newextgenes
    
    return(res)
    
  }
  
  extreadscounts1_adj <- extreadscounts1/(reads1depth/10^6)
  extreadscounts2_adj <- extreadscounts2/(reads2depth/10^6)
  
  extdepth1_adj <- sum(extreadscounts1_adj)
  extdepth2_adj <- sum(extreadscounts2_adj)
  
  weight1 <- extdepth1_adj/mean(c(extdepth1_adj, extdepth2_adj))
  weight2 <- extdepth2_adj/mean(c(extdepth1_adj, extdepth2_adj))
  
  extreadscounts1_adj <- extreadscounts1_adj/weight1
  extreadscounts2_adj <- extreadscounts2_adj/weight2
  
  #Add pseudo count!!!
  #extreadscounts1_adj <- extreadscounts1_adj + 10^-5/weight1
  #extreadscounts2_adj <- extreadscounts2_adj + 10^-5/weight2
  
  #Add pseudo count 1!!!
  extreadscounts1_adj <- extreadscounts1_adj + 1/weight1
  extreadscounts2_adj <- extreadscounts2_adj + 1/weight2
  
  extreadscounts_adj <- extreadscounts1_adj + extreadscounts2_adj
  
  
  
  extreadscounts_adj <- extreadscounts_adj/2
  
  scanextreads <- as.vector(extreadscounts_adj)
  
  windowsizes <- seq(from = 50, 
                     to = ceiling(length(scanextreads)/2/50)*50, 
                     by = 50)
  
  windowdiffs <- c()
  
  windowmeanses <- list()
  
  for(m in 1:length(windowsizes)){
    
    windowsize <- windowsizes[m]
    
    windowmeans1 <- c()
    windowmeans2 <- c()
    
    for(n in 1:length(scanextreads)){
      
      windowmean1 <- mean(scanextreads[n:min((windowsize + n - 1), 
                                             length(scanextreads))])
      windowmean2 <- mean(rev(scanextreads)[n:min((windowsize + n - 1), 
                                                  length(scanextreads))])
      
      windowmeans1 <- c(windowmeans1, windowmean1)
      windowmeans2 <- c(windowmeans2, windowmean2)
      
      
    }
    
    windowmeans <- (windowmeans1 + rev(windowmeans2))/2
    
    windowmeanses[[m]]  <- windowmeans
    
    windowdiff <- windowmeans[1] + windowmeans[length(windowmeans)] - 
      min(windowmeans)
    
    windowdiffs <- c(windowdiffs, windowdiff)
    
  }
  
  windowidxes <- order(-windowdiffs)
  
  defined <- FALSE
  
  for(windowidx in windowidxes){
    
    endpoint1 <- which.min(windowmeanses[[windowidx]]) - 1
    
    endpoint2 <- which.min(rev(windowmeanses[[windowidx]])) - 1
    
    if(endpoint1 <= ext1@ranges@width & endpoint2 <= ext2@ranges@width & 
       endpoint1 > 0 & endpoint2 > 0){
      
      defined <- TRUE
      
      break()
      
    }
    
  }
  
  if(defined == FALSE){
    
    windowidx <- windowidxes[1]
    
    endpoint1 <- which.min(windowmeanses[[windowidx]]) - 1
    
    endpoint2 <- which.min(rev(windowmeanses[[windowidx]])) - 1
    
    endpoint1 <- min(endpoint1, ext1@ranges@width)
    
    endpoint2 <- min(endpoint2, ext2@ranges@width)
    
  }
  
  endpoints <- c(endpoint1, endpoint2)
  
  res <- list()
  
  res$mixgenes <- c()
  
  if(sum(endpoints == 0) == 1){
    
    if(endpoints[1] == 0){
      
      return(NULL)
      
    }
    
    endpoint <- endpoints[endpoints != 0]
    strandsign <- strandsigns[endpoints != 0]
    ext <- exts[endpoints != 0]
    
    start <- ext@ranges@start
    end <- start + ext@ranges@width - 1
    
    if(strandsign == '+'){
      
      endpoint <- start + endpoint - 1
      
      endext <- GenomicRanges::GRanges(seqnames = chr, 
                                       ranges = IRanges::IRanges(start = start, 
                                                                 end = endpoint), 
                                       strand = strandsign)
      
      endext@elementMetadata <- ext@elementMetadata
      
    }else{
      
      endpoint <- end - endpoint + 1
      
      endext <- GenomicRanges::GRanges(seqnames = chr, 
                                       ranges = IRanges::IRanges(start = endpoint, 
                                                                 end = end), 
                                       strand = strandsign)
      
      endext@elementMetadata <- ext@elementMetadata
      
    }
    
    res$mixgenes <- endext
    
  }else{
    
    strandsign <- strandsigns[1]
    endpoint <- endpoints[1]
    ext <- exts[1]
    
    start <- ext@ranges@start
    end <- start + ext@ranges@width - 1
    
    if(strandsign == '+'){
      
      endpoint <- start + endpoint - 1
      
      endext <- GenomicRanges::GRanges(seqnames = chr, 
                                       ranges = IRanges::IRanges(start = start, 
                                                                 end = endpoint), 
                                       strand = strandsign)
      
      endext@elementMetadata <- ext@elementMetadata
      
    }else{
      
      endpoint <- end - endpoint + 1
      
      endext <- GenomicRanges::GRanges(seqnames = chr, 
                                       ranges = IRanges::IRanges(start = endpoint, 
                                                                 end = end), 
                                       strand = strandsign)
      
      endext@elementMetadata <- ext@elementMetadata
      
    }
    
    res$mixgenes <- c(res$mixgenes, endext)
    
    
  }
  
  if(length(res) == 0){
    
    res <- NULL
    
  }
  
  return(res)
  
}

smooth_fun <- function(x){
  
  datatype <- class(x)[1]
  
  if(datatype == 'Rle'){
    
    x <- as.vector(x)
    
  }
  
  res <- predict(
    
    locfit::locfit(x ~ locfit::lp(seq_along(x), nn = 0.1, h = 0.8)), 
    seq_along(x)
    
  )
  
  res[res < 0] <- 0
  
  if(datatype == 'Rle'){
    
    res <- S4Vectors::Rle(res)
    
  }
  
  return(res)
  
}

binratio <- function(i = 1, 
                     datalist,
                     window_num,
                     startshorten, 
                     endshorten,
                     reads1depth, 
                     reads2depth){
  
  gene <- datalist$genes[i]
  
  #return(gene)
  
  chr <- GenomicAlignments::seqlevelsInUse(gene)
  GenomeInfoDb::seqlevels(gene) <- chr
  #seqlevels contains all the chrs in a factor
  #seqlevelsInUse only contains the chrs present in x
  
  strandsign <- as.character(gene@strand@values)
  
  if(strandsign == '+'){
    
    tss <- gene@ranges@start
    tts <- tss + gene@ranges@width - 1
    
    inferranges <- GenomicRanges::GRanges(seqnames = chr, 
                                          ranges = IRanges::IRanges(start = tss, 
                                                                    end = tts), 
                                          strand = strandsign)
    
    genereads1 <- IRanges::subsetByOverlaps(x = datalist$subreads1, 
                                            ranges = inferranges)
    
    genereads2 <- IRanges::subsetByOverlaps(x = datalist$subreads2, 
                                            ranges = inferranges)
    
    GenomeInfoDb::seqlevels(genereads1) <- chr
    GenomeInfoDb::seqlevels(genereads2) <- chr
    
    fullreadscounts1 <- GenomicAlignments::coverage(genereads1, 
                                                    shift = -(tss - 1), 
                                                    width = (tts - tss + 1))
    #Calculate the coverage per base, like wig file, per chrom
    
    fullreadscounts2 <- GenomicAlignments::coverage(genereads2, 
                                                    shift = -(tss - 1), 
                                                    width = (tts - tss + 1))
    
    fullreadscounts1 <- fullreadscounts1[[1]]
    fullreadscounts2 <- fullreadscounts2[[1]]
    
    
    
    start <- gene@ranges@start + startshorten
    size <- floor((gene@ranges@width - startshorten - endshorten)/window_num)
    to <- start + size*window_num - 1
    #Important!
    #In R, the rule on range ends is different from HTSeq and most bed, gtf 
    #files. In these files, the start of a sequence is 0 and the end value on 
    #their record is actually not included in the range. 
    #However, in IRange and GRange, the sequence starts from 1 and the end value 
    #is included in the real sequence! So when import bed file into R 
    #session, it is necessary to adjust the coordinates to follow this R rule 
    #manually. R itself will not adjust it automatically.
    
    tss <- gene@ranges@start
    tts <- tss + gene@ranges@width - 1
    endshorten_adj <- tts - to
    
    genereadscounts1 <- fullreadscounts1[startshorten:(to - tss)]
    genereadscounts2 <- fullreadscounts2[startshorten:(to - tss)]
    
  }else{
    
    tss <- gene@ranges@start
    tts <- tss + gene@ranges@width - 1
    
    inferranges <- GenomicRanges::GRanges(seqnames = chr, 
                                          ranges = IRanges::IRanges(start = tss, 
                                                                    end = tts), 
                                          strand = strandsign)
    
    genereads1 <- IRanges::subsetByOverlaps(x = datalist$subreads1, 
                                            ranges = inferranges)
    
    genereads2 <- IRanges::subsetByOverlaps(x = datalist$subreads2, 
                                            ranges = inferranges)
    
    GenomeInfoDb::seqlevels(genereads1) <- chr
    GenomeInfoDb::seqlevels(genereads2) <- chr
    
    fullreadscounts1 <- GenomicAlignments::coverage(genereads1, 
                                                    shift = -(tss - 1), 
                                                    width = (tts - tss + 1))
    #Calculate the coverage per base, like wig file, per chrom
    
    fullreadscounts2 <- GenomicAlignments::coverage(genereads2, 
                                                    shift = -(tss - 1), 
                                                    width = (tts - tss + 1))
    
    fullreadscounts1[[1]] <- rev(fullreadscounts1[[1]])
    fullreadscounts2[[1]] <- rev(fullreadscounts2[[1]])
    
    fullreadscounts1 <- fullreadscounts1[[1]]
    fullreadscounts2 <- fullreadscounts2[[1]]
    
    
    
    to <- gene@ranges@start + gene@ranges@width - 1 - startshorten
    size <- floor((gene@ranges@width - startshorten - endshorten)/window_num)
    start <- to - size*window_num + 1
    #Important!
    #In R, the rule on range ends is different from HTSeq and most bed, gtf 
    #files. In these files, the start of a sequence is 0 and the end value on 
    #their record is actually not included in the range. 
    #However, in IRange and GRange, the sequence starts from 1 and the end value 
    #is included in the real sequence! So when import bed file into R 
    #session, it is necessary to adjust the coordinates to follow this R rule 
    #manually. R itself will not adjust it automatically.
    
    tss <- gene@ranges@start
    tts <- tss + gene@ranges@width - 1
    endshorten_adj <- start - tss
    
    genereadscounts1 <- fullreadscounts1[startshorten:(tts - start)]
    genereadscounts2 <- fullreadscounts2[startshorten:(tts - start)]
    
    
  }
  
  genedepth1 <- sum(genereadscounts1)
  genedepth2 <- sum(genereadscounts2)
  
  if(genedepth1 == 0 | genedepth2 == 0){
    return(NULL)
  }
  
  genereadscounts1_adj <- genereadscounts1/(reads1depth/10^6)
  genereadscounts2_adj <- genereadscounts2/(reads2depth/10^6)
  
  fullreadscounts1 <- fullreadscounts1/(reads1depth/10^6)
  fullreadscounts2 <- fullreadscounts2/(reads2depth/10^6)
  
  genedepth1_adj <- sum(genereadscounts1_adj)
  genedepth2_adj <- sum(genereadscounts2_adj)
  
  weight1 <- genedepth1_adj/mean(c(genedepth1_adj, genedepth2_adj))
  weight2 <- genedepth2_adj/mean(c(genedepth1_adj, genedepth2_adj))
  
  genereadscounts1_adj <- genereadscounts1_adj/weight1
  genereadscounts2_adj <- genereadscounts2_adj/weight2
  
  
  
  #Add pseudo count 1!!!
  #genereadscounts1_adj <- genereadscounts1_adj + 10^-5/weight1
  #genereadscounts2_adj <- genereadscounts2_adj + 10^-5/weight2
  
  #Add pseudo count!!!
  genereadscounts1_adj <- genereadscounts1_adj + 1/weight1
  genereadscounts2_adj <- genereadscounts2_adj + 1/weight2
  
  #Divide the genebody into 40 bins and calculate the read count of each bin
  
  starts <- seq(1, size*window_num, size)
  
  vi1 <- IRanges::Views(genereadscounts1_adj, start = starts, width = size)
  vi2 <- IRanges::Views(genereadscounts2_adj, start = starts, width = size)
  
  windows1 <- IRanges::viewSums(vi1)
  #viewSums calculates sums of the views in a Views or ViewsList object.
  #Here the total read num of each window can be calculated
  
  windows2 <- IRanges::viewSums(vi2)
  #viewSums calculates sums of the views in a Views or ViewsList object. 
  #Here the total read num of each window can be calculated
  
  #Generate the bin ratios between time points
  gene_ratio_adj <- windows2/windows1
  
  
  
  #gene_ratio_adj <- smooth_fun(x = gene_ratio_adj)
  
  
  
  if(sum(is.finite(gene_ratio_adj)) <= length(gene_ratio_adj)/2){
    
    return(NULL)
    
  }
  
  
  
  gene_ratio_adj <- data.frame(ratio_adj = gene_ratio_adj, 
                               stringsAsFactors = FALSE)
  
  res <- list(ratio_adj = gene_ratio_adj, 
              windowsize = size, 
              emis1_adj = genereadscounts1_adj, 
              emis2_adj = genereadscounts2_adj,
              emis1_complete = fullreadscounts1, 
              emis2_complete = fullreadscounts2, 
              endshorten_adj = endshorten_adj, 
              weight1 = weight1, 
              weight2 = weight2)
  
  return(res)
  
}



expandratio <- function(expand1 = binres$emis1_adj[startcoord:endcoord],
                        expand2 = binres$emis2_adj[startcoord:endcoord], 
                        
                        
                        
                        weight1 = weight1, 
                        weight2 = weight2){
  
  expand_ratio <- expand2/expand1
  
  #expand_ratio <- smooth_fun(x = expand_ratio)
  
  if(sum(is.finite(expand_ratio)) < length(expand_ratio)){
    
    #Add pseudo count!!!
    #expand_ratio <- (expand2 + 10^-5/weight2)/(expand1 + 10^-5/weight1)
    
    #Add pseudo count!!!
    expand_ratio <- (expand2 + 1/weight2)/(expand1 + 1/weight1)
    
  }
  
  expand_ratio <- data.frame(ratio_adj = expand_ratio, stringsAsFactors = FALSE)
  
  return(expand_ratio)
  
}

transpoint <- function(ratios, 
                       method = 'LSS', 
                       pythonpath = NULL, 
                       hmmseed = 2023){
  
  if(method == 'LSS'){
    
    expand_ratio <- ratios$ratio_adj
    
    startcoord <- 1
    endcoord <- length(expand_ratio)
    square_list <- c()
    
    for(k in startcoord:(endcoord - 1)){
      leftidx <- startcoord:k
      rightidx <- (k+1):endcoord
      #rightidx start from k+1, so the final transition point obtained from this function belong to the 1st state
      leftratio <- expand_ratio[leftidx - startcoord + 1]
      rightratio <- expand_ratio[rightidx - startcoord + 1]
      
      leftratio_clean <- leftratio#[(!is.na(leftratio)) & is.finite(leftratio)]
      rightratio_clean <- rightratio#[(!is.na(rightratio)) & is.finite(rightratio)]
      
      leftratio_var <- sum((leftratio_clean - mean(leftratio_clean))^2)
      rightratio_var <- sum((rightratio_clean - mean(rightratio_clean))^2)
      square <- leftratio_var + rightratio_var
      
      square_list <- c(square_list, square)
      
      
    }
    
    names(square_list) <- 1:length(square_list)
    minsquare <- min(square_list[(!is.na(square_list)) & is.finite(square_list)])
    
    if(is.infinite(minsquare)){
      
      return(NULL)
      
    }else{
      
      final_point <- as.numeric(names(square_list)[match(minsquare, as.vector(square_list))])
      
      fronts <- expand_ratio[1:final_point]
      latters <- expand_ratio[(final_point + 1):length(expand_ratio)]
      
      wilp <- wilcox.test(fronts, latters)$p.value
      
      frontmean <- mean(fronts)
      lattermean <- mean(latters)
      
      res <- list(point = final_point, front = frontmean, latter = lattermean, wilp = wilp)
      
      return(res)
      
    }
    
  }else{
      
    hmmpyfile <- system.file("python", "hmm_r.py", package = "proRate")
    #hmmpyfile <- 
    #  'C:/Users/Yu Liu/Desktop/proRate/proRate-master/inst/python/hmm_r.py'
    
    
    
    if(!is.null(pythonpath)){
      
      Sys.setenv(RETICULATE_PYTHON = pythonpath)
      
    }
    
    #reticulate::use_python(pydir)
    
    reticulate::py_config()
    
    reticulate::source_python(hmmpyfile)
    
    hmmres <- hmm_r(ratios = ratios, 
                    hmmseed = hmmseed)
    
    if(is.null(hmmres)){
      
      return(NULL)
      
    }
    
    final_point <- hmmres$point
    
    fronts <- hmmres$fronts
    
    latters <- hmmres$latters
    
    wilp <- wilcox.test(fronts, latters)$p.value
    
    frontmean <- mean(fronts)
    lattermean <- mean(latters)
    
    res <- list(point = final_point, front = frontmean, latter = lattermean, wilp = wilp)
    
    return(res)
    
    
  }
  
}

calsig <- function(emis1_adj = binres$emis1_adj, emis2_adj = binres$emis2_adj, point = distance){
  
  ratio_adj <- emis2_adj/emis1_adj
  
  activepart <- ratio_adj[1:point]
  inhibitpart <- ratio_adj[(point + 1):length(ratio_adj)]
  
  activepart <- activepart[is.finite(activepart)]
  inhibitpart <- inhibitpart[is.finite(inhibitpart)]
  
  wilp <- tryCatch({
    wilcox.test(activepart, inhibitpart)$p.value
  }, error = function(err){
    NA
  })
  
  activemean <- mean(activepart)
  inhibitmean <- mean(inhibitpart)
  
  res <- list(activeratiomean = activemean, inhibitratiomean = inhibitmean, wilp = wilp)
  return(res)
}



#reads1 = binres$ratio_adj
#reads2 = NULL
#endtype = endtypestr
#distance = point_idx$point
#genesym = gene$gene_id
#genelen = nrow(binres$ratio_adj)
#timediff = time
#saveplot = TRUE

#textsize = textsize
#titlesize = titlesize
#face = face

plotpairs <- function(reads1 = as.numeric(genefpm1), 
                      reads2 = as.numeric(genefpm2),
                      endtype = c('TSS', 'TTS'), 
                      distance, 
                      genesym = geneinfo$gene_id, 
                      genelen = geneinfo$width, 
                      timediff = time,
                      saveplot = FALSE, 
                      
                      textsize = 13,
                      titlesize = 15,
                      face = 'bold'){
  
  if(is.vector(reads1)){
    reads1 <- data.frame(reads = reads1, stringsAsFactors = FALSE)
  }
  
  row.names(reads1) <- 1:nrow(reads1)
  
  if(!is.null(reads2)){
    
    if(is.vector(reads2)){
      reads2 <- data.frame(reads = reads2, stringsAsFactors = FALSE)
    }
    row.names(reads2) <- 1:nrow(reads2)
  }
  
  start <- 1
  end <- nrow(reads1)
  distance <- distance
  
  points <- c(start, distance, end)
  pointlabels <- c(endtype[1], 'Transition', endtype[2])
  
  if(!is.null(time)){
    
    reads1$condition <- 'time2/time1'
    
  }else{
    
    reads1$condition <- 'condition/control'
    
  }
  
  reads1$xcoord <- as.numeric(row.names(reads1))
  reads_totals <- reads1
  ylabel <- 'Ratio'
  
  if(!is.null(reads2)){
    
    if(!is.null(time)){
      
      reads1$condition <- 'time1'
      reads2$condition <- 'time2'
      reads2$xcoord <- as.numeric(row.names(reads2))
      reads_totals <- rbind(reads1, reads2)
      ylabel <- 'FPM'
      reads_totals$condition <- factor(reads_totals$condition, levels = c('time2', 'time1'),
                                       ordered = TRUE)
      
    }else{
      
      reads1$condition <- 'control'
      reads2$condition <- 'condition'
      reads2$xcoord <- as.numeric(row.names(reads2))
      reads_totals <- rbind(reads1, reads2)
      ylabel <- 'FPM'
      reads_totals$condition <- factor(reads_totals$condition, levels = c('condition', 'control'),
                                       ordered = TRUE)
      
    }
    
    
  }
  
  names(reads_totals)[1] <- 'reads'
    
  if(!is.null(timediff)){
    
    plottitle <- paste0('(Distance = ', distance, ', Length = ', genelen, 
                        ', Time difference = ', timediff, 'min)')
    
  }else{
    
    plottitle <- paste0('(Distance = ', distance, ', Length = ', genelen, ')')
    
  }
  
  plottitle <- paste0(genesym, ' ', plottitle)
  
  #library(ggplot2)
  #library(scales)
  #library(grid)
  
  p <- ggplot2::ggplot(reads_totals, mapping = ggplot2::aes(x = xcoord, y = reads))
  
  if(!is.null(reads2)){
    p <- p + ggplot2::geom_line(size = 1, color = scales::hue_pal()(1))
  }else{
    p <- p + ggplot2::geom_point(size = 2, color = scales::hue_pal()(1))
  }
  
  p <- p + ggplot2::scale_x_continuous(breaks = points,
                                       labels = pointlabels) +
    ggplot2::xlab('') + ggplot2::ylab(ylabel) +
    ggplot2::geom_vline(xintercept = points[2], linetype = 2, color = 'red', size = 1) +
    ggplot2::facet_grid(condition~., scales = 'free_y') +
    ggplot2::ggtitle(plottitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust=1, angle = 90)) +
    ggplot2::theme(panel.spacing = grid::unit(0, 'lines')) + 
    
    
    
    ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face),
                   #plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                   axis.title = ggplot2::element_text(size = titlesize, face = face),
                   axis.text = ggplot2::element_text(size = textsize, face = face),
                   strip.text = ggplot2::element_text(size = textsize, face = face))
    
  
  if(saveplot == FALSE){
    print(p)
  }
  
  if(saveplot == TRUE){
    return(p)
  }
  
  
  
}



inferfunction <- function(i = 1, 
                          datalist, 
                          window_num, 
                          startshorten, 
                          endshorten, 
                          reads1depth, 
                          reads2depth, 
                          time, 
                          
                          method = 'LSS', 
                          pythonpath = NULL, 
                          hmmseed = 2023, 
                          
                          textsize = 13,
                          titlesize = 15,
                          face = 'bold'){
  
  binres <- binratio(i = i, 
                     datalist = datalist, 
                     window_num = window_num, 
                     startshorten = startshorten, 
                     endshorten = endshorten, 
                     reads1depth = reads1depth, 
                     reads2depth = reads2depth)
  
  weight1 <- binres$weight1
  weight2 <- binres$weight2
  
  if(is.null(binres)){
    #next()
    return(NULL)
  }
  
  endshorten <- binres$endshorten_adj
  
  
  
  ###################
  #Start 1st infer
  
  #LSS method
  point_idx <- transpoint(ratios = binres$ratio_adj, 
                          method = method, 
                          pythonpath = pythonpath, 
                          hmmseed = hmmseed)
  
  
  
  if(is.null(point_idx)){
    #next()
    return(NULL)
  }
  
  gene <- datalist$genes[i]
  
  if(!is.null(time)){
    
    endtypestr <- c(paste0('TSS+', startshorten, 'bp'), 
                    paste0('TTS-', endshorten, 'bp'))
    
  }else{
    
    endtypestr <- c(paste0('ExonStart+', startshorten, 'bp'), 
                    paste0('DistPolyA-', endshorten, 'bp'))
    
  }
  
  binresplot <- plotpairs(reads1 = binres$ratio_adj, 
                          reads2 = NULL,
                          endtype = endtypestr,
                          distance = point_idx$point, 
                          genesym = gene$gene_id,
                          genelen = nrow(binres$ratio_adj), 
                          timediff = time, 
                          saveplot = TRUE, 
                          
                          textsize = textsize,
                          titlesize = titlesize,
                          face = face)
  
  ###################
  forward <- round(binres$windowsize * (point_idx$point - 1))
  backward <- round(binres$windowsize * (point_idx$point + 1))
  
  #focus on the single window with transition point
  
  startcoord <- 1 + forward
  endcoord <- backward
  
  expand_ratio <- expandratio(expand1 = binres$emis1_adj[startcoord:endcoord],
                              expand2 = binres$emis2_adj[startcoord:endcoord], 
                              weight1 = weight1, 
                              weight2 = weight2)
  
  
  
  if(is.null(expand_ratio)){
    
    final_point <- list(point = binres$windowsize, front = NA, latter = NA, wilp = NA)
    
  }else{
    
    
    
    #Start 2nd infer
    
    #LSS method
    final_point <- transpoint(ratios = expand_ratio, 
                              method = method, 
                              pythonpath = pythonpath, 
                              hmmseed = hmmseed)
    
    if(is.null(final_point)){
      
      final_point <- list(point = binres$windowsize, front = NA, latter = NA, wilp = NA)
      
    }
    
    
    
  }
  
  if(!is.null(time)){
    
    endtypestr <- c(paste0('TSS+', forward + startshorten, 'bp'), 
                    paste0('TSS+', (backward - 1) + startshorten, 'bp'))
    
  }else{
    
    endtypestr <- c(paste0('ExonStart+', forward + startshorten, 'bp'), 
                    paste0('ExonStart+', (backward - 1) + startshorten, 'bp'))
    
  }
  
  expandplot <- plotpairs(reads1 = expand_ratio$ratio_adj, 
                          reads2 = NULL,
                          endtype = endtypestr,
                          distance = final_point$point, 
                          genesym = gene$gene_id,
                          genelen = length(expand_ratio$ratio_adj), 
                          timediff = time,
                          saveplot = TRUE, 
                          
                          textsize = textsize,
                          titlesize = titlesize,
                          face = face)
  
  res <- list()
  
  res$gene_id <- gene$gene_id
  
  res$distance <- startshorten + forward + final_point$point
  
  res$wilp1 <- point_idx$wilp
  res$wilp2 <- final_point$wilp
  
  res$activemean1 <- point_idx$front
  res$inhibitmean1 <- point_idx$latter
  
  res$activemean2 <- final_point$front
  res$inhibitmean2 <- final_point$latter
  
  res$singlegeneread1 <- binres$emis1_complete
  res$singlegeneread2 <- binres$emis2_complete
  res$binresplot <- binresplot
  res$expandplot <- expandplot
  
  
  subreport <- data.frame(gene_id = res$gene_id,
                          distance = res$distance,
                          frontbinratio = res$activemean1,
                          latterbinratio = res$inhibitmean1,
                          diffbinratio = res$inhibitmean1 - res$activemean1,
                          binpval = res$wilp1,
                          frontextendratio = res$activemean2,
                          latterextendratio = res$inhibitmean2,
                          diffextendratio = res$inhibitmean2 - res$activemean2,
                          extendpval = res$wilp2,
                          stringsAsFactors = FALSE)
  
  reslist <- list(subreport = subreport,
                  singlegenereads1 = res$singlegeneread1,
                  singlegenereads2 = res$singlegeneread2,
                  binresplotlist = res$binresplot,
                  expandplotlist = res$expandplot)
  
  
  return(reslist)
  
  
}

#time1file = wt0file
#time2file = wt15file
#targetfile = NULL
#genomename = 'mm10'
#time = 15
#strandmethod = 1
#threads = 1

#lencutoff = 40000
#fpkmcutoff = 1

#startshorten = 1000
#endshorten = 1000
#window_num = 40
#method = 'LSS'
#pythonpath = 'C:/Users/yuabr/anaconda3/python.exe'
#difftype = 1
#utr = FALSE
#utrexts = NULL

#textsize = 13
#titlesize = 15
#face = 'bold'

#'Calculate gene elongation rate from a pair of Pro-seq or Gro-seq data
#'
#'Calculate gene elongation rate from a pair of Pro-seq or Gro-seq data, using 
#'the LSS (least sum of squares) or HMM (hidden Markov model) method. 
#'
#'@param time1file The reference Pro-seq/Gro-seq bam file, corresponding to 
#'  the experimental condition of no transcriptional inhibitor treatment. Can 
#'  be a string indicating the directory of the file, or a GAlignmentPairs 
#'  object, a GAlignments object, or a GRanges object from the original bam 
#'  file.
#'@param time2file The treatment Pro-seq/Gro-seq bam file, corresponding to 
#'  the experimental condition of transcriptional inhibitor treatment for a 
#'  specific time (e.g. DRB treatment for 15 min). Can be a string indicating 
#'  the directory of the bam file, or a GAlignmentPairs object, a GAlignments 
#'  object, or a GRanges object from the original file.
#'@param targetfile A txt file with the genes whose transcriptional rates need 
#'  to be calculated. Should contain columns named as chr, start, end, strand, 
#'  and gene_id. It can also be NULL, so that the genes in the genome set by 
#'  the parameter \code{genomename} will be analyzed. However, in any case, 
#'  the genes should have a length longer than the one set by the parameter 
#'  \code{lencutoff}, and also longer than the one of 2*(\code{startshorten} + 
#'  \code{endshorten}), which is set by the parameters \code{startshorten} and 
#'  \code{endshorten}.
#'@param gene_ids A vector with gene symbols indicating the ones need to be 
#'  analyzed. In addition to \code{targetefile} and \code{genomename}, this 
#'  parameter also indicates the genes to be analyzed. The final ones should 
#'  belong to the intersection of these parameters, and they also need to have 
#'  a length longer than the one set by the parameter \code{lencutoff}, and 
#'  also longer than the one of 2*(\code{startshorten} + \code{endshorten}), 
#'  which is set by the parameters \code{startshorten} and \code{endshorten}. 
#'  Can also be NULL, so that no restriction will be added from it. 
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetfile} is NULL. Can be "mm10" for mouse or "hg38" for 
#'  human.
#'@param time An integer indicating the inhibitor treatment time difference 
#'  between the \code{time1file} and the \code{time2file}, using min as its 
#'  unit.
#'@param strandmethod Indicate the strand specific method used when preparing 
#'  the sequencing library, can be 1 for the directional ligation method, 2 
#'  for the dUTP method, and 0 for a non-strand specific library. In addition, 
#'  if the sample is sequenced using a single strand method, set it as 3.
#'@param threads Number of threads to perform parallelization. Default is 1.
#'@param lencutoff The cutoff on gene length (bp). Only genes longer than this 
#'  cutoff can be considered for analysis. Default is 70000.
#'@param fpkmcutoff The cutoff value on gene FPKM. Only genes with an FPKM 
#'  value greater than the cutoff in \code{time1file} can be considered for 
#'  analysis. Default is 1.
#'@param startshorten Before inferring a gene's transcription rate, its first 
#'  1000 bp (or other length) and last 1000 bp (or other length) regions will 
#'  be discarded to avoid the unstable reads at the transcription starting and 
#'  ending stages. However, these regions' lengths can be changed by setting 
#'  this parameter \code{startshorten} and the other \code{endshorten}. This 
#'  one is used to set the length of the transcription starting region. Its 
#'  default value is 1000, so that the first 1000 bp region will be discarded.
#'@param endshorten Before inferring a gene's transcription rate, its first 
#'  1000 bp (or other length) and last 1000 bp (or other length) regions will 
#'  be discarded to avoid the unstable reads at the transcription starting and 
#'  ending stages. However, these regions' lengths can be changed by setting 
#'  this parameter \code{endshorten} and the other \code{startshorten}. This 
#'  one is used to set the length of the transcription ending region. Default 
#'  is 1000, so that the last 1000 bp region will be discarded.
#'@param window_num Before inferring a gene's transcription rate, the function 
#'  will divide this gene into 40 bins (or other bin number). For each bin, 
#'  the normalized read count ratio between the treatment and the reference 
#'  files will be calculated, so a vector with 40 ratios (or other bin number) 
#'  will be generated. Then, the LSS or HMM method will be used to find the 
#'  transition bin between the gene's transcription inhibited region and the 
#'  normal reads region. After that, this identified transition bin and its 
#'  downstream neighbor will be merged and expanded to the single-base-pair 
#'  level, and the LSS or HMM method will be further used on them to find the 
#'  transition base pair in this region. The parameter \code{window_num} here 
#'  is used to set the bin number to be divided for each gene. Default value 
#'  is 40.
#'@param method The method to be used for transcription rate inference. The 
#'  default value is "LSS", so that the least sum of squares method will be 
#'  used. Can also be "HMM", so that the hidden Markov model will be used. 
#'@param pythonpath The HMM method is base on \code{Python}, so the directory 
#'  of the \code{Python} interpreter you want to use should be transferred to 
#'  the function via this parameter, and two \code{Python} modules should be 
#'  installed in your \code{Python} environment, including \code{numpy} and 
#'  \code{hmmlearn}.
#'@param hmmseed The HMM method involves random processes, so a random seed 
#'  should be set via this parameter to repeat the results. Default value is 
#'  1234, can also be other integers, such as 2023.
#'@param difftype In most cases, the treatment and reference Pro-seq/Gro-seq 
#'  files are from experiments treating cells with transcription inhibitors, 
#'  such as DRB (5,6-dichloro-1-beta-d-ribofuranosylbenzimidazole), so that 
#'  the normal transcription will be repressed for a specific time, generating 
#'  a reads-depleted region upstream of the normal transcription region. For 
#'  such inhibitor-based experiments, this parameter \code{difftype} should be 
#'  set as 1. However, in some cases, the treatment and reference Pro-seq/Gro- 
#'  seq files can also come from experiments treating cells with transcription 
#'  activators, e.g., treating MCF-7 human breast cancer cells with E2 (17- 
#'  beta-estradiol), making the reads-depleted region downstream, rather than 
#'  upstream, of the normal transcription region, which is in contrast to the 
#'  DRB (inhibitor) experiments. For such activator-based experiments, this 
#'  parameter should be set as 2. In addition, this function \code{calrate} 
#'  can also infer proximal polyA alternative sites for genes, and in this 
#'  case, the parameter \code{time1file} should be an RNA-seq file with genes 
#'  using distal polyA sites; the parameter \code{time2file} needs an RNA-seq 
#'  file with genes using proximal polyA sites; another parameter \code{utr} 
#'  should be set as TRUE; and the parameter \code{difftype} here should be 
#'  set as 2. The default value of \code{difftype} is 1.
#'@param utr In addition to inferring transcription rates from Pro-seq/Gro-seq 
#'  data, \code{calrate} can also infer proximal polyA alternative sites for 
#'  genes. In this case, the parameter \code{time1file} should be an RNA-seq 
#'  file with genes using distal polyA sites; the parameter \code{time2file} 
#'  needs an RNA-seq file with genes using proximal polyA sites; the parameter 
#'  \code{difftype} should be set as 2; and the current parameter \code{utr} 
#'  should be set as TRUE. The default value of \code{utr} is FALSE, so the 
#'  function will perform inference on transcription rates, not on proximal 
#'  polyA sites.
#'@param utrexts When the former parameter \code{utr} is set as TRUE to infer 
#'  proximal polyA sites for genes, this parameter can be used to provide a 
#'  txt file with the genes' last exons whose proximal polyA sites need to be 
#'  identified. Should contain columns named as chr, start, end, strand, and 
#'  gene_id. It can also be set as NULL, so that the genes' last exons in the 
#'  genome set by the parameter \code{genomename} will be analyzed. However, 
#'  in the latter case, the original last exons from the genome will be first 
#'  adjusted so that for the ones with a length > 10000 bp, the proximal polyA 
#'  sites will be inferred directly in them, but for the ones with a length <= 
#'  10000 bp, their lengths will be extended to 10000 bp first, and then the 
#'  proximal sites will be identified within the extended exons. On the other 
#'  hand, if the last exons are provided with this parameter \code{utrexts}, 
#'  they will never be extended to 10000 bp, and the proximal polyA inference 
#'  will be performed directly on them. In addition to the proximal sites, the 
#'  distal ones will also be defined by the function from the \code{time1file} 
#'  and \code{time2file}'s reads. It is performed with a sliding window method 
#'  on the last exon regions, and this step is before the proximal polyA sites 
#'  inference but after the exon extension step. It should also be noted that 
#'  the last exons should have an FPKM value in the \code{time1file} greater 
#'  than the cutoff set by \code{fpkmcutoff}.
#'@param textsize In addition to returning a data frame to show the inference 
#'  results, this function will also generate several plots to show them, and 
#'  the font size for the plot texts is set by this parameter. Default is 13.
#'@param titlesize The font size for the plot titles. Default is 15.
#'@param face The font face for the plot texts. Default is "bold".
#'@return A list including a slot named "report", which is a data frame with 
#'  the inferred transcription elongation rates, or genes' proximal and distal 
#'  polyA sites, as well as other information, such as the genes' coordinates, 
#'  the results' significance, etc. In addition, the result list also contains 
#'  other slots, such as "binplots" and "expandplots", which contains the data 
#'  that can be used to plot the inference results.
#'
#'
#'
#'@examples
#'library(proRate)
#'
#'wt0file <- system.file("extdata", "wt0.bam", package = "proRate")
#'wt15file <- system.file("extdata", "wt15.bam", package = "proRate")
#'
#'wtrates <- calrate(time1file = wt0file, 
#'                   time2file = wt15file, 
#'                   time = 15, 
#'                   strandmethod = 1, 
#'                   
#'                   genomename = "mm10", 
#'                   lencutoff = 40000, 
#'                   fpkmcutoff = 1, 
#'                   
#'                   threads = 1, 
#'                   
#'                   startshorten = 1000, 
#'                   endshorten = 1000, 
#'                   window_num = 40, 
#'                   
#'                   method = "LSS", 
#'                   pythonpath = NULL, 
#'                   
#'                   difftype = 1)
#'
#'
#'
#'@export
calrate <- function(time1file, 
                    time2file,
                    targetfile = NULL, 
                    gene_ids = NULL, 
                    genomename = 'mm10',
                    time,
                    strandmethod = 1, 
                    threads = 1,
                    
                    lencutoff = 70000,
                    fpkmcutoff = 1, 
                    
                    startshorten = 1000, 
                    endshorten = 1000, 
                    window_num = 40, 
                    method = 'LSS', 
                    pythonpath = NULL, 
                    hmmseed = 1234,  
                    #hmmseed = 2023,  
                    difftype = 1, 
                    utr = FALSE, 
                    utrexts = NULL, 
                    
                    textsize = 13,
                    titlesize = 15,
                    face = 'bold'){
  
  if(strandmethod %in% c(0, 1, 2)){
    singleend <- FALSE
  }else{
    singleend <- TRUE
    strandmethod <- 4
  }
  
  if(utr == TRUE){
    
    startshorten <- 0
    endshorten <- 0
    difftype <- 2
    time <- NULL
    
  }
  
  ##########
  #Libraries required
  #library(GenomicAlignments)
  
  if(!is.null(targetfile) & utr == FALSE){
      
    genes <- read.table(targetfile, sep = '\t', header = TRUE,
                        stringsAsFactors = FALSE, quote = '',
                        check.names = FALSE)
    
    genes <- genes[(abs(genes$end - genes$start) + 1) >= 
                     max(2*(startshorten + endshorten), lencutoff), , 
                   drop = FALSE]
    
    if(nrow(genes) == 0){
      
      return(NULL)
      
    }
    
    genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
    names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
    
    genes <- GenomicRanges::GRanges(seqnames = genes$seqnames,
                                    ranges = IRanges::IRanges(start = genes$start, end = genes$end),
                                    strand = genes$strand,
                                    gene_id = genes$gene_id)
    genes <- GenomicAlignments::sort(genes)
    
  }else{
    
    genes <- get(paste0(genomename, '.packagegenes'))
    #genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    genes <- GenomicRanges::GRanges(genes)
    
    if(utr == FALSE){
      
      genes <- genes[genes@ranges@width >= max(2*(startshorten + endshorten), lencutoff)]
      
    }
    
    if(length(genes) == 0){
      
      return(NULL)
      
    }
    
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- GenomeInfoDb::keepSeqlevels(x = genes, 
                                         value = unique(as.character(genes@seqnames)),
                                         pruning.mode = 'coarse')
    #Keeps only the seqlevels in `value` and removes all others.
    #x Any object having a Seqinfo class in which the seqlevels will be kept.
    #value A named or unnamed character vector.
    
    if(utr == FALSE){
      
      genes$updis <- NULL
      genes$dndis <- NULL
      
    }
    
  }
  
  if(!is.null(targetfile) & utr == FALSE){
    
    genesdis <- GenomicRanges::distanceToNearest(genes, ignore.strand = TRUE)
    
    mixgenes <- genesdis[genesdis@elementMetadata$distance == 0]
    
    mixgenecoords <- NULL
    
    if(length(mixgenes) > 0){
      
      mixgenes1 <- genes$gene_id[mixgenes@from]
      mixgenes2 <- genes$gene_id[mixgenes@to]
      
      mixgenes1strands <- as.vector(genes@strand)[mixgenes@from]
      mixgenes2strands <- as.vector(genes@strand)[mixgenes@to]
      
      mixgenes <- data.frame(mixgenes1 = mixgenes1, 
                             mixgenes2 = mixgenes2, 
                             stringsAsFactors = FALSE)
      
      mixgenes <- mixgenes[mixgenes1strands != mixgenes2strands, , drop = FALSE]
      
      if(nrow(mixgenes) > 0){
        
        row.names(mixgenes) <- 1:nrow(mixgenes)
        
      }
      
    }
    
    genecoords <- genes$gene_id
    
    if(sum(nrow(mixgenes) > 0)){
      
      mixgenecoords <- mixgenes
      mixgenes <- mixgenecoords
      
      mixgenecoords <- unique(c(mixgenecoords$mixgenes1, mixgenecoords$mixgenes2))
      
      genecoords <- setdiff(genecoords, mixgenecoords)
      
    }
    
    genecoords <- genes[genes$gene_id %in% genecoords]
    
    if(length(genecoords) == 0){
      
      return(NULL)
      
    }else{
      
      genecoords <- GenomicAlignments::sort(genecoords)
      
      names(genecoords) <- 1:length(genecoords)
      
      genes <- genecoords
      
    }
    
    
  }
  
  if(!is.null(gene_ids)){
    
    genes <- genes[genes$gene_id %in% gene_ids]
    
    names(genes) <- 1:length(genes)
    
  }
  
  ###################
  
  #Read in the strand-specific paired-end data
  
  reads1 <- parsereadsfile(readsfile = time1file, strandmethod = strandmethod)
  reads2 <- parsereadsfile(readsfile = time2file, strandmethod = strandmethod)
  
  reads1depth <- length(reads1)
  reads2depth <- length(reads2)
  
  
  #genes <- genes[genes$gene_id %in% 
  #                 c('Cpt1a', 'Kmt5b', 'Fxn', 'Csf2ra', 'Tesmin', 'Gal')]
  #names(genes) <- 1:length(genes)
  
  
  if(utr == TRUE){
    
    threads <- max(1, round(threads))
    
    #Infer distal polyA site
    if(is.null(utrexts)){
      
      #lastexons <- makeutrlist(genes = genes, genomename = genomename)
      
      lastexons <- get(paste0('lastexons.', genomename))
      
      #lastexons <- readRDS(paste0('lastexons.', genomename, '.rds'))
      
      nullgenes <- setdiff(genes$gene_id, lastexons$gene_id)
      
      nullgenes <- genes[genes$gene_id %in% nullgenes]
      
      genes <- genes[!(genes$gene_id %in% nullgenes$gene_id)]
      
      lastexons <- lastexons[match(genes$gene_id, lastexons$gene_id)]
      
      trimeddndis <- 10000 - lastexons@ranges@width
      trimeddndis[trimeddndis < 0] <- 0
      trimeddndis[trimeddndis > genes$dndis] <- genes$dndis[trimeddndis > genes$dndis]
      
      
      pluslimits <- genes@ranges@start + genes@ranges@width - 1 + trimeddndis
      minuslimits <- genes@ranges@start - trimeddndis
      
      plusends <- pluslimits
      minusends <- minuslimits
      
      extends <- plusends
      extends[as.vector(lastexons@strand) == '-'] <- 
        minusends[as.vector(lastexons@strand) == '-']
      
      plusstarts <- lastexons@ranges@start
      minusstarts <- lastexons@ranges@start + lastexons@ranges@width - 1
      
      extstarts <- plusstarts
      extstarts[as.vector(lastexons@strand) == '-'] <- 
        minusstarts[as.vector(lastexons@strand) == '-']
      
      genes$extstarts <- extstarts
      genes$extends <- extends
      
      genes$extstarts[as.vector(genes@strand) == '-'] <- 
        extends[as.vector(genes@strand) == '-']
      genes$extends[as.vector(genes@strand == '-')] <- 
        extstarts[as.vector(genes@strand) == '-']
      
      exts <- GenomicRanges::GRanges(seqnames = genes@seqnames, 
                                     ranges = IRanges::IRanges(start = genes$extstarts, 
                                                               end = genes$extends), 
                                     strand = genes@strand, 
                                     gene_id = genes$gene_id)
      
      names(exts) <- 1:length(exts)
      
      extsdis <- GenomicRanges::distanceToNearest(exts, ignore.strand = TRUE)
      
      mixexts <- extsdis[extsdis@elementMetadata$distance == 0]
      
      mixgenes <- NULL
      
      if(length(mixexts) > 0){
        
        mixexts1 <- exts$gene_id[mixexts@from]
        mixexts2 <- exts$gene_id[mixexts@to]
        
        mixexts1strands <- as.vector(exts@strand)[mixexts@from]
        mixexts2strands <- as.vector(exts@strand)[mixexts@to]
        
        mixexts <- data.frame(mixexts1 = mixexts1, 
                              mixexts2 = mixexts2, 
                              stringsAsFactors = FALSE)
        
        mixexts <- mixexts[mixexts1strands != mixexts2strands, , drop = FALSE]
        
        if(nrow(mixexts) > 0){
          
          row.names(mixexts) <- 1:nrow(mixexts)
          
        }
        
      }
      
      if(!is.null(targetfile)){
          
        genes <- read.table(targetfile, sep = '\t', header = TRUE,
                            stringsAsFactors = FALSE, quote = '',
                            check.names = FALSE)
        
        genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
        names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
        
        extgenes <- exts$gene_id[exts$gene_id %in% genes$gene_id]
        
        if(sum(nrow(mixexts) > 0)){
          
          mixgenes <- mixexts[mixexts$mixexts1 %in% genes$gene_id, , drop = FALSE]
          mixexts <- mixgenes
          
          mixgenes <- unique(c(mixgenes$mixexts1, mixgenes$mixexts2))
          
          extgenes <- setdiff(extgenes, mixgenes)
          
        }
        
        
        
      }else{
        
        extgenes <- exts$gene_id
        
        if(sum(nrow(mixexts) > 0)){
          
          mixgenes <- mixexts
          mixexts <- mixgenes
          
          mixgenes <- unique(c(mixgenes$mixexts1, mixgenes$mixexts2))
          
          extgenes <- setdiff(extgenes, mixgenes)
          
        }
        
      }
      
      extgenes <- exts[exts$gene_id %in% extgenes]
      
      mixgenes <- exts[exts$gene_id %in% mixgenes]
      
      endexts <- GenomicRanges::GRanges(seqnames = character(),
                                        ranges = IRanges::IRanges(start = numeric(), 
                                                                  end = numeric()),
                                        strand = character(),
                                        gene_id = character())
      
      if(length(mixgenes) > 0){
        
        mixgenes <- IRanges::subsetByOverlaps(x = mixgenes, 
                                              ranges = range(c(reads1, reads2)))
        
      }
      
      if(length(mixgenes) > 0){
        
        mixgenes <- GenomicAlignments::sort(mixgenes)
        
        names(mixgenes) <- 1:length(mixgenes)
        
        nseq <- seq(1, length(mixgenes), 1)
        
        if(threads == 1){
          
          mixreslist <- list()
          n <- 1
          for(n in nseq){
            
            mixres <- mixscan(i = n, 
                              mixext = mixgenes, 
                              extpair = mixexts, 
                              
                              reads1 = reads1, 
                              reads2 = reads2, 
                              reads1depth = reads1depth, 
                              reads2depth = reads2depth)
            
            mixreslist[[n]] <- mixres
            
            cat(paste0('n = ', n, '\n'))
            
          }
          
          
        }else{
          
          #library(doParallel)
          #threads <- 8
          
          cores <- parallel::detectCores()
          cl <- parallel::makeCluster(min(min(threads, length(mixgenes)), cores))
          
          doParallel::registerDoParallel(cl)
          
          date()
          `%dopar%` <- foreach::`%dopar%`
          
          mixreslist <- foreach::foreach(i = nseq, 
                                         .export = ls(name = globalenv()), .packages = c('GenomicRanges')) %dopar% {
                                           #.export = NULL) %dopar% {
                                           mixscan(i = i, 
                                                   mixext = mixgenes, 
                                                   extpair = mixexts, 
                                                   
                                                   reads1 = reads1, 
                                                   reads2 = reads2, 
                                                   reads1depth = reads1depth, 
                                                   reads2depth = reads2depth)
                                         }
          date()
          
          parallel::stopCluster(cl)
          
          unregister_dopar()
          
          
        }
        
        for(n in 1:length(mixreslist)){
          
          if(!is.null(names(mixreslist[[n]]))){
            
            if(names(mixreslist[[n]]) == 'extgenes'){
              
              extgenes <- c(extgenes, mixreslist[[n]]$extgenes)
              
            }else{
              
              endexts <- c(endexts, mixreslist[[n]]$mixgenes)
              
            }
            
          }
          
        }
        
        if(length(extgenes) > 0){
          extgenes <- unique(extgenes)
          names(extgenes) <- 1:length(extgenes)
        }
        
        if(length(endexts) > 0){
          endexts <- unique(endexts)
          names(endexts) <- 1:length(endexts)
        }
        
      }
      
      if(length(extgenes) > 0){
        
        extgenes <- IRanges::subsetByOverlaps(x = extgenes, 
                                              ranges = range(c(reads1, reads2)))
        
      }
      
      if(length(extgenes) > 0){
        
        extgenes <- GenomicAlignments::sort(extgenes)
        
        names(extgenes) <- 1:length(extgenes)
        
        nseq <- seq(1, length(extgenes), 1)
        
        if(threads == 1){
          
          endextlist <- list()
          n <- 1
          for(n in nseq){
            
            endext <- dirscan(i = n, 
                              ext = extgenes, 
                              reads1 = reads1, 
                              reads2 = reads2, 
                              reads1depth = reads1depth, 
                              reads2depth = reads2depth)
            
            endextlist[[n]] <- endext
            
            cat(paste0('n = ', n, '\n'))
            
          }
          
          
        }else{
          
          #library(doParallel)
          #threads <- 8
          
          cores <- parallel::detectCores()
          cl <- parallel::makeCluster(min(min(threads, length(extgenes)), cores))
          
          doParallel::registerDoParallel(cl)
          
          date()
          `%dopar%` <- foreach::`%dopar%`
          
          endextlist <- foreach::foreach(i = nseq, 
                                         .export = ls(name = globalenv()), .packages = c('GenomicRanges')) %dopar% {
                                           #.export = NULL) %dopar% {
                                           dirscan(i = i, 
                                                   ext = extgenes, 
                                                   reads1 = reads1, 
                                                   reads2 = reads2, 
                                                   reads1depth = reads1depth, 
                                                   reads2depth = reads2depth)
                                           }
          date()
          
          parallel::stopCluster(cl)
          
          unregister_dopar()
          
          
        }
        
        endextlist <- Filter(Negate(is.null), endextlist)
        
        endextlist <- unlist(endextlist)
        
        if(is.list(endextlist)){
          endextlist <- do.call(c, endextlist)
          
          if(length(endextlist) > 0){
            names(endextlist) <- 1:length(endextlist)
          }
          
        }else{
          endextlist <- NULL
        }
        
        
        endexts <- c(endexts, endextlist)
        
      }
      
      if(length(endexts) == 0){
        
        return(NULL)
        
      }else{
        
        endexts <- endexts[endexts$gene_id %in% exts$gene_id]
        
        if(length(endexts) == 0){
          
          return(NULL)
          
        }
        
      }
      
      exts <- endexts
      
      exts <- GenomicAlignments::sort(exts)
      
      names(exts) <- 1:length(exts)
      
    }else{
        
      exts <- read.table(utrexts, sep = '\t', header = TRUE,
                         stringsAsFactors = FALSE, quote = '',
                         check.names = FALSE)
      
      exts <- exts[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(exts) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      exts <- exts[exts$gene_id %in% genes$gene_id]
      
      exts <- GenomicRanges::GRanges(seqnames = exts$seqnames,
                                     ranges = IRanges::IRanges(start = exts$start, 
                                                               end = exts$end),
                                     strand = exts$strand,
                                     gene_id = exts$gene_id)
      exts <- GenomicAlignments::sort(exts)
      
      names(exts) <- 1:length(exts)
      
      extsdis <- GenomicRanges::distanceToNearest(exts, ignore.strand = TRUE)
      
      mixexts <- extsdis[extsdis@elementMetadata$distance == 0]
      
      mixgenes <- NULL
      
      if(length(mixexts) > 0){
        
        mixexts1 <- exts$gene_id[mixexts@from]
        mixexts2 <- exts$gene_id[mixexts@to]
        
        mixexts1strands <- as.vector(exts@strand)[mixexts@from]
        mixexts2strands <- as.vector(exts@strand)[mixexts@to]
        
        mixexts <- data.frame(mixexts1 = mixexts1, 
                              mixexts2 = mixexts2, 
                              stringsAsFactors = FALSE)
        
        mixexts <- mixexts[mixexts1strands != mixexts2strands, , drop = FALSE]
        
        if(nrow(mixexts) > 0){
          
          row.names(mixexts) <- 1:nrow(mixexts)
          
        }
        
      }
      
      extgenes <- exts$gene_id
      
      if(sum(nrow(mixexts) > 0)){
        
        mixgenes <- mixexts
        mixexts <- mixgenes
        
        mixgenes <- unique(c(mixgenes$mixexts1, mixgenes$mixexts2))
        
        extgenes <- setdiff(extgenes, mixgenes)
        
      }
      
      extgenes <- exts[exts$gene_id %in% extgenes]
      
      mixgenes <- exts[exts$gene_id %in% mixgenes]
      
      endexts <- GenomicRanges::GRanges(seqnames = character(),
                                        ranges = IRanges::IRanges(start = numeric(), 
                                                                  end = numeric()),
                                        strand = character(),
                                        gene_id = character())
      
      
      if(length(mixgenes) > 0){
        
        mixgenes <- IRanges::subsetByOverlaps(x = mixgenes, 
                                              ranges = range(c(reads1, reads2)))
        
      }
      
      if(length(mixgenes) > 0){
        
        mixgenes <- GenomicAlignments::sort(mixgenes)
        
        names(mixgenes) <- 1:length(mixgenes)
        
        nseq <- seq(1, length(mixgenes), 1)
        
        if(threads == 1){
          
          mixreslist <- list()
          n <- 1
          for(n in nseq){
            
            mixres <- mixscan(i = n, 
                              mixext = mixgenes, 
                              extpair = mixexts, 
                              
                              reads1 = reads1, 
                              reads2 = reads2, 
                              reads1depth = reads1depth, 
                              reads2depth = reads2depth)
            
            mixreslist[[n]] <- mixres
            
            cat(paste0('n = ', n, '\n'))
            
          }
          
          
        }else{
          
          #library(doParallel)
          #threads <- 8
          
          cores <- parallel::detectCores()
          cl <- parallel::makeCluster(min(min(threads, length(mixgenes)), cores))
          
          doParallel::registerDoParallel(cl)
          
          date()
          `%dopar%` <- foreach::`%dopar%`
          
          mixreslist <- foreach::foreach(i = nseq, 
                                         .export = ls(name = globalenv()), .packages = c('GenomicRanges')) %dopar% {
                                           #.export = NULL) %dopar% {
                                           mixscan(i = i, 
                                                   mixext = mixgenes, 
                                                   extpair = mixexts, 
                                                   
                                                   reads1 = reads1, 
                                                   reads2 = reads2, 
                                                   reads1depth = reads1depth, 
                                                   reads2depth = reads2depth)
                                         }
          
          date()
          
          parallel::stopCluster(cl)
          
          unregister_dopar()
          
          
        }
        
        for(n in 1:length(mixreslist)){
          
          if(!is.null(names(mixreslist[[n]]))){
            
            if(names(mixreslist[[n]]) == 'extgenes'){
              
              extgenes <- c(extgenes, mixreslist[[n]]$extgenes)
              
            }else{
              
              endexts <- c(endexts, mixreslist[[n]]$mixgenes)
              
            }
            
          }
          
        }
        
        
        if(length(extgenes) > 0){
          extgenes <- unique(extgenes)
          names(extgenes) <- 1:length(extgenes)
        }
        
        if(length(endexts) > 0){
          endexts <- unique(endexts)
          names(endexts) <- 1:length(endexts)
        }
        
      }
      
      if(length(extgenes) > 0){
        
        extgenes <- IRanges::subsetByOverlaps(x = extgenes, 
                                              ranges = range(c(reads1, reads2)))
        
      }
      
      if(length(extgenes) > 0){
        
        extgenes <- GenomicAlignments::sort(extgenes)
        
        names(extgenes) <- 1:length(extgenes)
        
        endexts <- c(endexts, extgenes)
        
      }
      
      if(length(endexts) == 0){
        
        return(NULL)
        
      }else{
        
        endexts <- endexts[endexts$gene_id %in% exts$gene_id]
        
        if(length(endexts) == 0){
          
          return(NULL)
          
        }
        
      }
      
      exts <- endexts
      
      exts <- GenomicAlignments::sort(exts)
      
      names(exts) <- 1:length(exts)
      
    }
    
    #FPKM screen###
    exts <- getfpkm(features = exts, reads = reads1, readsdepth = reads1depth, 
                    singleend = singleend)
    exts <- exts[!is.na(exts$fpkm)]
    exts <- exts[exts$fpkm != 0]
    if(length(exts) == 0){
      return(NULL)
    }
    exts <- GenomicAlignments::sort(exts)
    names(exts) <- 1:length(exts)
    
    if(!is.null(fpkmcutoff)){
      
      exts <- exts[exts$fpkm > fpkmcutoff]
      if(length(exts) == 0){
        return(NULL)
      }
      exts <- GenomicAlignments::sort(exts)
      names(exts) <- 1:length(exts)
      
    }
    
    #Prepare for parallel######
    
    if(threads > length(exts)){
      threads <- length(exts)
    }
    
    
    extranges <- range(exts)
    
    extreads1 <- IRanges::subsetByOverlaps(x = reads1, ranges = extranges)
    
    extreads2 <- IRanges::subsetByOverlaps(x = reads2, ranges = extranges)
    
    extreads <- c(extreads1, extreads2)
    extreadsranges <- range(extreads)
    exts <- IRanges::subsetByOverlaps(x = exts, ranges = extreadsranges)
    
    exts <- GenomicAlignments::sort(exts)
    extreads1 <- GenomicAlignments::sort(extreads1)
    extreads2 <- GenomicAlignments::sort(extreads2)
    
    GenomeInfoDb::seqlevels(extreads1) <- 
      GenomicAlignments::seqlevelsInUse(extreads1)
    
    GenomeInfoDb::seqlevels(extreads2) <- 
      GenomicAlignments::seqlevelsInUse(extreads2)
    
    if(length(exts) == 0 | length(extreads1) == 0 | length(extreads2) == 0){
      return(NULL)
    }
    
    names(exts) <- 1:length(exts)
    
    datalist <- list(genes = exts, subreads1 = extreads1, subreads2 = extreads2)
    
    
    
  }else{
    
    #FPKM screen###
    genes <- getfpkm(features = genes, reads = reads1, readsdepth = reads1depth, 
                     singleend = singleend)
    genes <- genes[!is.na(genes$fpkm)]
    genes <- genes[genes$fpkm != 0]
    if(length(genes) == 0){
      return(NULL)
    }
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    if(!is.null(fpkmcutoff)){
      
      genes <- genes[genes$fpkm > fpkmcutoff]
      if(length(genes) == 0){
        return(NULL)
      }
      genes <- GenomicAlignments::sort(genes)
      names(genes) <- 1:length(genes)
      
    }
    
    #Prepare for parallel######
    
    threads <- max(1, round(threads))
    
    if(threads > length(genes)){
      threads <- length(genes)
    }
    
    generanges <- range(genes)
    
    subreads1 <- IRanges::subsetByOverlaps(x = reads1, ranges = generanges)
    #subsetByOverlaps returns the subset of `x` that has an overlap hit with a range 
    #in `ranges` using the specified `findOverlaps` parameters.
    
    subreads2 <- IRanges::subsetByOverlaps(x = reads2, ranges = generanges)
    #subsetByOverlaps returns the subset of `x` that has an overlap hit with a range 
    #in `ranges` using the specified `findOverlaps` parameters.
    
    subreads <- c(subreads1, subreads2)
    subreadsranges <- range(subreads)
    genes <- IRanges::subsetByOverlaps(x = genes, ranges = subreadsranges)
    #subsetByOverlaps returns the subset of `x` that has an overlap hit with a range 
    #in `ranges` using the specified `findOverlaps` parameters.
    
    if(length(genes) == 0){
      
      return(NULL)
      
    }
    
    rm(subreads)
    
    genes <- GenomicAlignments::sort(genes)
    subreads1 <- GenomicAlignments::sort(subreads1)
    subreads2 <- GenomicAlignments::sort(subreads2)
    
    GenomeInfoDb::seqlevels(subreads1) <- 
      GenomicAlignments::seqlevelsInUse(subreads1)
    #seqlevels contains all the chrs in a factor
    #seqlevelsInUse only contains the chrs present in x
    
    GenomeInfoDb::seqlevels(subreads2) <- 
      GenomicAlignments::seqlevelsInUse(subreads2)
    #seqlevels contains all the chrs in a factor
    #seqlevelsInUse only contains the chrs present in x
    
    if(length(genes) == 0 | length(subreads1) == 0 | length(subreads2) == 0){
      return(NULL)
    }
    
    names(genes) <- 1:length(genes)
    
    datalist <- list(genes = genes, subreads1 = subreads1, subreads2 = subreads2)
    
    
  }
  
  
  geneframe <- as.data.frame(datalist$genes)
  
  for(i in 1:ncol(geneframe)){
    
    if(is.factor(geneframe[,i])){
      geneframe[,i] <- as.character(geneframe[,i])
    }
    
  }
  
  row.names(geneframe) <- 1:nrow(geneframe)
  
  
  iseq <- seq(1, length(datalist$genes), 1)
  
  gc()
  
  if(threads == 1){
    
    reslists <- list()
    i <- 22
    for(i in iseq){
      
      reslist <- inferfunction(i = i, 
                               datalist = datalist, 
                               window_num = window_num, 
                               startshorten = startshorten, 
                               endshorten = endshorten, 
                               reads1depth = reads1depth, 
                               reads2depth = reads2depth, 
                               time = time, 
                               
                               method = method, 
                               pythonpath = pythonpath, 
                               hmmseed = hmmseed, 
                               
                               textsize = textsize,
                               titlesize = titlesize,
                               face = face)
      
      reslists[[i]] <- reslist
      
    }
    
    
  }else{
    
    #library(doParallel)
    #threads <- 8
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    doParallel::registerDoParallel(cl)
    
    date()
    `%dopar%` <- foreach::`%dopar%`
    reslists <- foreach::foreach(i = iseq, 
                                 .export = ls(name = globalenv()), .packages = c('GenomicRanges')) %dopar% {
                                 #.export = NULL, .packages = c('GenomicRanges')) %dopar% {
                                 inferfunction(i = i, 
                                               datalist = datalist, 
                                               window_num = window_num, 
                                               startshorten = startshorten, 
                                               endshorten = endshorten, 
                                               reads1depth = reads1depth, 
                                               reads2depth = reads2depth, 
                                               time = time, 
                                               
                                               method = method, 
                                               pythonpath = pythonpath, 
                                               hmmseed = hmmseed, 
                                               
                                               textsize = textsize,
                                               titlesize = titlesize,
                                               face = face)
                                   }
    #.export character vector of variables to export. This can be useful when 
    #accessing a variable that isn't defined in the current environment. The 
    #default valule is NULL
    
    #.packages charactor vector of packages that the tasks depend on. If the 
    #expression, ex, requires an R package to be loaded, this option can be used 
    #to load that package on each of the works.
    
    date()
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
    
  }
  
  report <- NULL
  genereads1 <- list()
  genereads2 <- list()
  binplots <- list()
  expandplots <- list()
  
  for(i in 1:length(reslists)){
    
      report <- rbind(report, reslists[[i]]$subreport)
      genereads1[[i]] <- reslists[[i]]$singlegenereads1
      genereads2[[i]] <- reslists[[i]]$singlegenereads2
      binplots[[i]] <- reslists[[i]]$binresplotlist
      expandplots[[i]] <- reslists[[i]]$expandplotlist
      
  }
  
  
  
  if(is.null(report) & length(genereads1) == 0 & length(genereads2) == 0 & 
     length(binplots) == 0 & length(expandplots) == 0){
    
    res <- list()
    res[['report']] <- report
    res[['genereads1']] <- genereads1
    res[['genereads2']] <- genereads2
    
    res[['binplots']] <- binplots
    res[['expandplots']] <- expandplots
    
    return(res)
  }
  
  
  
  genereads1 <- Filter(Negate(is.null), genereads1)
  genereads2 <- Filter(Negate(is.null), genereads2)
  binplots <- Filter(Negate(is.null), binplots)
  expandplots <- Filter(Negate(is.null), expandplots)
  
  names(genereads1) <- names(genereads2) <- 
    names(binplots) <- names(expandplots) <- report$gene_id
  
  
  if(utr == FALSE){
    
    report$time <- time
    report$rate <- report$distance/report$time
    
  }
  
  report$binpadj <- p.adjust(p = report$binpval, method = 'BH')
  report$extendpadj <- p.adjust(p = report$extendpval, method = 'BH')
  
  
  report$significance <- 'nonsignificant'
  
  if(difftype == 1){
    
    report$significance[report$diffbinratio > 0 & report$binpadj < 0.05] <- 'significant'
    
  }else if(difftype == 2){
    
    report$significance[report$diffbinratio < 0 & report$binpadj < 0.05] <- 'significant'
    
  }else{
    
    report$significance[report$binpadj < 0.05] <- 'significant'
    
  }
  
  
  geneframe <- geneframe[c('gene_id', 'seqnames', 'start', 'end', 'strand',
                           setdiff(names(geneframe), c('gene_id', 'seqnames', 'start', 'end', 'strand')))]
  
  names(report)[1] <- 'gene_id'
  
  
  sigreport <- subset(report, significance == 'significant')
  nonsigreport <- subset(report, significance == 'nonsignificant')
  
  sigreport <- sigreport[order(sigreport$binpadj, sigreport$binpval,
                               sigreport$extendpadj, sigreport$extendpval,
                               -abs(sigreport$diffbinratio), -abs(sigreport$diffextendratio),
                               -sigreport$distance),]
  nonsigreport <- nonsigreport[order(nonsigreport$binpadj, nonsigreport$binpval,
                                     nonsigreport$extendpadj, nonsigreport$extendpval,
                                     -abs(nonsigreport$diffextendratio), -abs(nonsigreport$diffextendratio),
                                     -nonsigreport$distance),]
  report <- rbind(sigreport, nonsigreport)
  row.names(report) <- 1:nrow(report)
  
  geneframe <- subset(geneframe, gene_id %in% report$gene_id)
  
  geneframe <- geneframe[match(report$gene_id, geneframe$gene_id),]
  row.names(geneframe) <- 1:nrow(geneframe)
  report <- cbind(geneframe, report[-1])
  
  names(report)[2] <- 'chr'
  
  if(utr == FALSE){
    
    maincols <- c('gene_id',
                  'distance', 'time', 'rate',
                  'significance', 'binpadj', 'binpval',
                  'frontbinratio', 'latterbinratio', 'diffbinratio',
                  'chr', 'start', 'end', 'strand',
                  'extendpadj', 'extendpval',
                  'frontextendratio', 'latterextendratio', 'diffextendratio')
    
  }else{
    
    report$ProxPolyA <- report$start + report$distance - 1
    report$ProxPolyA[report$strand == '-'] <- 
      (report$end - report$distance + 1)[report$strand == '-']
    report$ProxPolyA <- paste0(report$chr, ':', report$ProxPolyA)
    
    report$DistPolyA <- report$end
    report$DistPolyA[report$strand == '-'] <- report$start[report$strand == '-']
    report$DistPolyA <- paste0(report$chr, ':', report$DistPolyA)
    
    maincols <- c('gene_id',
                  'distance', 'ProxPolyA', 'DistPolyA',
                  'significance', 'binpadj', 'binpval',
                  'frontbinratio', 'latterbinratio', 'diffbinratio',
                  'chr', 'start', 'end', 'strand',
                  'extendpadj', 'extendpval',
                  'frontextendratio', 'latterextendratio', 'diffextendratio')
    
  }
  
  othercols <- setdiff(names(report), maincols)
  allcols <- c(maincols, othercols)
  report <- report[allcols]
  row.names(report) <- 1:nrow(report)
  
  names(report)[names(report) == 'width'] <- 'genewidth'
  
  if(utr == TRUE){
    
    names(report)[names(report) == 'start'] <- 'ExtStart'
    names(report)[names(report) == 'end'] <- 'ExtEnd'
    names(report)[names(report) == 'genewidth'] <- 'ExtWidth'
    
  }
  
  res <- list()
  res[['report']] <- report
  res[['genereads1']] <- genereads1
  res[['genereads2']] <- genereads2
  
  res[['binplots']] <- binplots
  res[['expandplots']] <- expandplots
  
  return(res)
  
}



#lssres <- calrate(time1file = wt0file, 
#                  time2file = wt15file, 
#                  time = 15, 
#                  strandmethod = 1, 
#                  targetfile = NULL, 
#                  genomename = 'mm10', 
#                  lencutoff = 40000, 
#                  fpkmcutoff = 1, 
#                
#                  threads = 1, 
#                  
#                  startshorten = 1000, 
#                  endshorten = 1000, 
#                  window_num = 40, 
#                  
#                  method = 'LSS', 
#                  pythonpath = 'C:/Users/yuabr/anaconda3/python.exe', 
#                  difftype = 1)



#lssres <- calrate(time1file = fp0file, 
#                  time2file = fp12file, 
#                  time = 12, 
#                  strandmethod = 4, 
#                  targetfile = NULL, 
#                  genomename = 'mm10', 
#                  lencutoff = 40000, 
#                  fpkmcutoff = 1, 
#                  
#                  threads = 1, 
#                  
#                  startshorten = 1000, 
#                  endshorten = 1000, 
#                  window_num = 40, 
#                  
#                  method = 'LSS', 
#                  pythonpath = 'C:/Users/yuabr/anaconda3/python.exe', 
#                  difftype = 1, 
#                  utrexts = NULL, 
#                  
#                  textsize = 13, 
#                  titlesize = 15, 
#                  face = 'bold')



#hmmres <- calrate(time1file = fp0file, 
#                  time2file = fp12file, 
#                  time = 12, 
#                  strandmethod = 4, 
#                  targetfile = NULL, 
#                  genomename = 'mm10', 
#                  lencutoff = 40000, 
#                  fpkmcutoff = 1, 
#                  
#                  threads = 1, 
#                  
#                  startshorten = 1000, 
#                  endshorten = 1000, 
#                  window_num = 40, 
#                  
#                  method = 'HMM', 
#                  pythonpath = 'C:/Users/yuabr/anaconda3/python.exe', 
#                  difftype = 1, 
#                  utrexts = NULL, 
#                  
#                  textsize = 13, 
#                  titlesize = 15, 
#                  face = 'bold')





#'Calculate gene elongation rate for multiple pairs of Pro-seq or Gro-seq data 
#'
#'Calculate gene elongation rate for multiple pairs of Pro-seq or Gro-seq data 
#'with the LSS (least sum of squares) or HMM (hidden Markov model) method. 
#'
#'@param time1files The reference Pro-seq/Gro-seq bam files, corresponding to 
#'  the experimental condition of no transcriptional inhibitor treatment. Can 
#'  be a vector with elements as strings indicating the directories of the bam 
#'  files.
#'@param time2files The treatment Pro-seq/Gro-seq bam files, corresponding to 
#'  the treatment of transcriptional inhibitor for specific times (e.g. DRB 
#'  treatment for 15 min, 30 min, etc). Should be a vector with elements as 
#'  strings indicating the directories of the bam files.
#'@param targetfile A txt file with the genes whose transcriptional rates need 
#'  to be calculated. Should contain columns named as chr, start, end, strand, 
#'  and gene_id. It can also be NULL, so that the genes in the genome set by 
#'  the parameter \code{genomename} will be analyzed. However, in any case, 
#'  the genes should have a length longer than the one set by the parameter 
#'  \code{lencutoff}, and also longer than the one of 2*(\code{startshorten} + 
#'  \code{endshorten}), which is set by the parameters \code{startshorten} and 
#'  \code{endshorten}.
#'@param gene_ids A vector with gene symbols indicating the ones need to be 
#'  analyzed. In addition to \code{targetefile} and \code{genomename}, this 
#'  parameter also indicates the genes to be analyzed. The final ones should 
#'  belong to the intersection of these parameters, and they also need to have 
#'  a length longer than the one set by the parameter \code{lencutoff}, and 
#'  also longer than the one of 2*(\code{startshorten} + \code{endshorten}), 
#'  which is set by the parameters \code{startshorten} and \code{endshorten}. 
#'  Can also be NULL, so that no restriction will be added from it. 
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetfile} is NULL. Can be "mm10" for mouse or "hg38" for 
#'  human.
#'@param times The treatment time differences between the \code{time1files} 
#'  and the \code{time2files}, using min as their units. Should be a vector 
#'  with each element as the time difference between the matched elements in 
#'  the \code{time1files} vector and the \code{time2files} vector.
#'@param strandmethod Indicate the strand specific method used when preparing 
#'  the sequencing libraries, can be 1 for the directional ligation method, 2 
#'  for the dUTP method, and 0 for non-strand specific librares. In addition, 
#'  if the samples are sequenced using a single strand method, set it as 3.
#'@param threads Number of threads to do the parallelization. Default is 1.
#'@param mergerefs Whether to merge all the reference data contained in the 
#'  \code{time1files} vector to one, and then use it as a unified reference 
#'  for all the \code{time2files}. Default is TRUE.
#'@param mergecases Whether to merge all the treatment data contained in the 
#'  \code{time2files} vector to one. Default is FALSE.
#'@param lencutoff The cutoff on gene length (bp). Only genes longer than this 
#'  cutoff can be considered for analysis. Default is 70000.
#'@param fpkmcutoff The cutoff value on gene FPKM. Only genes with an FPKM 
#'  value greater than the cutoff in the reference data can be considered for 
#'  analysis. Default is 1.
#'@param startshorten Before inferring a gene's transcription rate, its first 
#'  1000 bp (or other length) and last 1000 bp (or other length) regions will 
#'  be discarded to avoid the unstable reads at the transcription starting and 
#'  ending stages. However, these regions' lengths can be changed by setting 
#'  this parameter \code{startshorten} and the other \code{endshorten}. This 
#'  one is used to set the length of the transcription starting region. Its 
#'  default value is 1000, so that the first 1000 bp region will be discarded.
#'@param endshorten Before inferring a gene's transcription rate, its first 
#'  1000 bp (or other length) and last 1000 bp (or other length) regions will 
#'  be discarded to avoid the unstable reads at the transcription starting and 
#'  ending stages. However, these regions' lengths can be changed by setting 
#'  this parameter \code{endshorten} and the other \code{startshorten}. This 
#'  one is used to set the length of the transcription ending region. Default 
#'  is 1000, so that the last 1000 bp region will be discarded.
#'@param window_num Before inferring a gene's transcription rate, the function 
#'  will divide this gene into 40 bins (or other bin number). For each bin, 
#'  the normalized read count ratio between the treatment and the reference 
#'  files will be calculated, so a vector with 40 ratios (or other bin number) 
#'  will be generated. Then, the LSS or HMM method will be used to find the 
#'  transition bin between the gene's transcription inhibited region and the 
#'  normal reads region. After that, this identified transition bin and its 
#'  downstream neighbor will be merged and expanded to the single-base-pair 
#'  level, and the LSS or HMM method will be further used on them to find the 
#'  transition base pair in this region. The parameter \code{window_num} here 
#'  is used to set the bin number to be divided for each gene. Default value 
#'  is 40.
#'@param method The method to be used for transcription rate inference. The 
#'  default value is "LSS", so that the least sum of squares method will be 
#'  used. Can also be "HMM", so that the hidden Markov model will be used. 
#'@param pythonpath The HMM method is base on \code{Python}, so the directory 
#'  of the \code{Python} interpreter you want to use should be transferred to 
#'  the function via this parameter, and two \code{Python} modules should be 
#'  installed to your \code{Python} environment, including \code{numpy} and 
#'  \code{hmmlearn}.
#'@param hmmseed The HMM method involves random processes, so a random seed 
#'  should be set via this parameter to repeat the results. Default value is 
#'  1234, can also be other integers, such as 2023.
#'@param difftype In most cases, the treatment and reference Pro-seq/Gro-seq 
#'  files are from experiments treating cells with transcription inhibitors, 
#'  such as DRB (5,6-dichloro-1-beta-d-ribofuranosylbenzimidazole), so that 
#'  the normal transcription will be repressed for a specific time, generating 
#'  a reads-depleted region upstream of the normal transcription region. For 
#'  such inhibitor-based experiments, this parameter \code{difftype} should be 
#'  set as 1. However, in some cases, the treatment and reference Pro-seq/Gro- 
#'  seq files can also come from experiments treating cells with transcription 
#'  activators, e.g., treating MCF-7 human breast cancer cells with E2 (17- 
#'  beta-estradiol), making the reads-depleted region downstream, rather than 
#'  upstream, of the normal transcription region, which is in contrast to the 
#'  DRB (inhibitor) experiments. For such activator-based experiments, this 
#'  parameter should be set as 2. In addition, this function \code{mcalrate} 
#'  can also infer proximal polyA alternative sites for genes, and to perform 
#'  this analysis, the parameter \code{time1files} needs to contain RNA-seq 
#'  files with genes using distal polyA sites; the parameter \code{time2files} 
#'  should have RNA-seq files with genes using proximal polyA sites; another 
#'  parameter \code{utr} needs to be TRUE; and the parameter \code{difftype} 
#'  here should be set as 2. The default value of \code{difftype} is 1.
#'@param utr In addition to inferring transcription rates from Pro-seq/Gro-seq 
#'  data, \code{mcalrate} can also infer proximal polyA alternative sites for 
#'  genes. In this case, the parameter \code{time1files} should be an RNA-seq 
#'  files with genes using distal polyA sites; the parameter \code{time2files} 
#'  needs RNA-seq files with genes using proximal polyA sites; the parameter 
#'  \code{difftype} should be set as 2; and the current parameter \code{utr} 
#'  should be set as TRUE. The default value of \code{utr} is FALSE, so the 
#'  function will perform inference on transcription rates, not on proximal 
#'  polyA sites.
#'@param utrexts When the former parameter \code{utr} is set as TRUE to infer 
#'  proximal polyA sites for genes, this parameter can be used to provide a 
#'  txt file with the genes' last exons whose proximal polyA sites need to be 
#'  identified. Should contain columns named as chr, start, end, strand, and 
#'  gene_id. It can also be set as NULL, so that the genes' last exons in the 
#'  genome set by the parameter \code{genomename} will be analyzed. However, 
#'  in the latter case, the original last exons from the genome will be first 
#'  adjusted so that for the ones with a length > 10000 bp, the proximal polyA 
#'  sites will be inferred directly in them, but for the ones with a length <= 
#'  10000 bp, their lengths will be extended to 10000 bp first, and then the 
#'  proximal sites will be identified within the extended exons. On the other 
#'  hand, if the last exons are provided with this parameter \code{utrexts}, 
#'  they will never be extended to 10000 bp, and the proximal polyA inference 
#'  will be performed directly on them. In addition to the proximal sites, the 
#'  distal polyA sites will also be defined by the function from the Pro-seq/
#'  Gro-seq pairs defined by \code{time1files} and \code{time2files}. It is 
#'  performed with a sliding window method on the last exons, and it is before 
#'  the proximal polyA sites inference but after the exon extension step. It 
#'  should be noted that for a specific Pro-seq/Gro-seq file included in the 
#'  parameter \code{time1files}, it also set a cutoff for the last exons to be 
#'  analyzed with the polyA sites inference, i.e., their FPKM values in this 
#'  file should be greater than \code{fpkmcutoff}. 
#'@param textsize In addition to returning a data frame to show the inference 
#'  results, this function will also generate several plots to show them, and 
#'  the font size for the plot texts is set by this parameter. Default is 13.
#'@param titlesize The font size for the plot titles. Default is 15.
#'@param face The font face for the plot texts. Default is "bold".
#'@return A list with several sub-lists and each of them includes a slot named 
#'  "report", which is a data frame with the inferred transcription rates, or 
#'  genes' proximal and distal polyA sites, as well as other information, such 
#'  as the genes' coordinates, the results' significance, etc. A sub-list also 
#'  contains other slots, such as "binplots" and "expandplots", which contains 
#'  the data that can be used to plot the inference results.
#'@export
mcalrate <- function(time1files,
                     time2files,
                     targetfile = NULL, 
                     gene_ids = NULL, 
                     genomename = 'mm10',
                     times,
                     strandmethod = 0, 
                     threads = 1,
                     mergerefs = TRUE, 
                     mergecases = FALSE,
                     
                     lencutoff = 70000,
                     fpkmcutoff = 1, 
                     
                     startshorten = 1000, 
                     endshorten = 1000, 
                     window_num = 40, 
                     method = 'LSS', 
                     pythonpath = NULL, 
                     hmmseed = 1234,  
                     #hmmseed = 2023,  
                     difftype = 1, 
                     
                     utr = FALSE, 
                     utrexts = NULL, 
                     
                     textsize = 13,
                     titlesize = 15,
                     face = 'bold'){
  
  time1files <- getreadslist(timefiles = time1files, mergefiles = mergerefs, strandmode = strandmethod)
  time2files <- getreadslist(timefiles = time2files, mergefiles = mergecases, strandmode = strandmethod)
  
  timefilelength <- length(time2files)
  
  reslist <- list()
  i <- 1
  for(i in 1:timefilelength){
    
    time2file <- time2files[[i]]
    time1file <- time1files[[min(length(time1files), i)]]
    
    if(utr == TRUE){
      
      time <- NULL
      
    }else{
      
      time <- times[i]
      
    }
    
    res <- calrate(time1file = time1file, 
                   time2file = time2file,
                   targetfile = targetfile, 
                   gene_ids = gene_ids, 
                   genomename = genomename,
                   time = time,
                   strandmethod = strandmethod, 
                   threads = threads,
                   
                   lencutoff = lencutoff,
                   fpkmcutoff = fpkmcutoff, 
                   
                   startshorten = startshorten, 
                   endshorten = endshorten, 
                   window_num = window_num, 
                   method = method, 
                   pythonpath = pythonpath, 
                   hmmseed = hmmseed, 
                   difftype = difftype, 
                   
                   utr = utr, 
                   utrexts = utrexts, 
                   
                   textsize = textsize,
                   titlesize = titlesize,
                   face = face)
    
    reslist[[i]] <- res
    
  }
  
  return(reslist)
  
}


