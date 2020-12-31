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



#'Generate a GRanges object of UCSC known genes in sepcific genome
#'
#'Generate a GRanges object of UCSC known genes in specific genome.
#'
#'@param genomename The name of the genome, can be 'mm10' for mouse or 'hg38' 
#'  for human.
#'@param save Whether save the result GRanges as a file in work directory, 
#'  default is TRUE.
#'@return A GRanges object of UCSC known genes, with several gene information 
#'  provided, including seqnames, ranges, strands, gene symbols, distance to 
#'  the nearest upstream gene, distance to the nearest downstream gene, gene 
#'  GC content, and gene exon coverage. Genes overlapping with each other will 
#'  not be included.
#'@export
makegenelist <- function(genomename = 'mm10', save = TRUE){
  
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
  
  
  genesyms <- genesyms[complete.cases(genesyms),]
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
  
  genecoords <- genegccont(origenes = genecoords, genomename = genomename)
  genecoords <- geneexonfreq(origenes = genecoords, genomename = genomename)
  
  genecoords <- GenomicAlignments::sort(genecoords)
  genes <- genecoords
  names(genes) <- 1:length(genes)
  
  if(save == TRUE){
    saveRDS(genes, paste0(genomename, '.packagegenes.rds'))
  }
  
  return(genes)
  
}



getfpkm <- function(features = genes, reads = reads1,
                    readsdepth = reads1depth){
  
  fcounts <- GenomicAlignments::summarizeOverlaps(features = features, reads = reads,
                                                  mode = 'IntersectionNotEmpty',
                                                  singleEnd = FALSE, ignore.strand = FALSE)
  #Add singleEnd = FALSE for pair-end
  #Add ignore.strand = FALSE for strand specific
  
  fcounts <- as.data.frame(SummarizedExperiment::assays(fcounts)$counts)
  featurewidth <- features@ranges@width
  
  factorM <- readsdepth/10^6
  factorK <- featurewidth/10^3
  
  fpkms <- fcounts/factorK/factorM
  features$fpkm <- fpkms$reads
  
  return(features)
  
}

parsereadsfile <- function(readsfile = time1file, strandmethod = 1){
  
  #library(GenomicAlignments)
  
  if(is.character(readsfile)){
    if(file.exists(readsfile)){
      
      if(strandmethod != 0){
        reads_ori <- GenomicAlignments::readGAlignmentPairs(readsfile, strandMode = strandmethod)
      }else{
        reads_ori <- GenomicAlignments::readGAlignments(readsfile)
      }
      reads <- as(reads_ori, 'GRanges')
      
    }
  }else if(class(readsfile)[1] == 'GAlignmentPairs'){
    
    reads <- as(readsfile, 'GRanges')
    
  }else if(class(readsfile)[1] == 'GRanges'){
    
    reads <- readsfile
    
  }
  
  return(reads)
  
}

getreadslist <- function(timefiles = time1files, mergefiles = mergerefs, strandmode = strandmethod){
  
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

towig <- function(reads, strand = "*", chrom = NULL){
  
  if(!is.null(chrom)){
    reads <- reads[GenomicAlignments::seqnames(reads) == chrom,]
  }
  
  
  readsList <- S4Vectors::split(reads, GenomicAlignments::seqnames(reads))
  #important!
  #split divides the data in a vector-like object 'x' into the
  #groups defined by 'f'.
  #Here split the bam file data into different chroms
  
  #Change reads strand
  readsList <- S4Vectors::endoapply(readsList, function(x){
    
    if(strand == "*"){
      BiocGenerics::strand(x) <- "*"
    }else{
      x <- x[BiocGenerics::strand(x) == strand,]
    }
    
    return(x)
  })
  
  #extract the specific strand reads defined by the fuction parameter
  
  #lapply returns a list of the same length as X, each element of which is
  #the result of applying FUN to the corresponding element of X.
  
  #sapply is a user-friendly version and wrapper of lapply by default
  #returning a vector, matrix.
  
  #endoapply performs the endomorphic equivalents of lapply
  #by returning objects of the same class as the inputs rather than a list.
  
  testlist <- as.list(readsList)
  
  #Select reads for specific normal chrs
  if(!is.null(chrom)){
    testlist <- testlist[chrom]
  }else{
    chroms <- names(testlist)
    normchroms <- chroms[grepl('chr', chroms)]
    normchroms <- unique(normchroms)
    testlist <- testlist[normchroms]
    normchroms <- names(testlist)[lapply(testlist, length) > 0]
    testlist <- testlist[normchroms]
  }
  
  i <- 1
  H <- list()
  for(i in 1:length(testlist)){
    x <- testlist[[i]]
    GenomeInfoDb::seqlevels(x) <- GenomicAlignments::seqlevelsInUse(x)
    #seqlevels contains all the chrs in a factor
    #seqlevelsInUse only contains the chrs present in x
    
    cov <- GenomicAlignments::coverage(x)[[1]]
    #Important
    #Caluculate the coverage per base, like wig file, per chrom
    to <- length(cov)
    #length(cov) = seqlengths(x), it is chrom length
    starts <- 1:to
    vi <- XVector::Views(cov, start=starts, width=1)
    
    H[[i]] <- S4Vectors::Rle(XVector::viewSums(vi))
    #viewSums calculate sums of the views in a Views or ViewsList object.
    #Here the total read num of each window can be calculated
    #Rle(values, lengths): This constructor creates an Rle instance out of
    #an atomic vector or factor object 'values' and an integer or
    #numeric vector 'lengths' with all positive elements that
    #represent how many times each value is repeated. The length
    #of these two vectors must be the same.
    
    names(H)[i] <- names(testlist)[i]
  }
  
  
  return(H)
  
}

binratio <- function(singlegeneinfo,
                     Fp1 = NULL, Fp2 = NULL, Fm1 = NULL, Fm2 = NULL, window_num = window_num,
                     startshorten = startshorten, endshorten = endshorten,
                     depth1 = depth1, depth2 = depth2, savegenenames = savegenesyms){
  
  if(singlegeneinfo[1, 5] == '+'){
    
    start <- singlegeneinfo[1, 2] + startshorten
    end   <- singlegeneinfo[1, 3] - endshorten
    
    #Important!
    #In R, the rule on range ends is different from HTSeq and most bed, gtf
    #files. In these files, the start of a sequence is 0 and the end value on
    #their record is actually not included in the range.
    #However, in IRange and GRange, the sequence start from 1 and the end value
    #is included in the real sequence! So when import bed file into R
    #session, it is necessary to adjust the coordiante to follow this R rule
    #manully. R itself will not adjust it automatically.
    
    #Extract the read count of the range from start to end, each number of
    #the returned value represents the read count in each bp
    emis1  <-
      (as.numeric(Fp1[[as.character(singlegeneinfo[1, 1])]]))[c(start:end)]
    emis2  <-
      (as.numeric(Fp2[[as.character(singlegeneinfo[1, 1])]]))[c(start:end)]
    
    if(singlegeneinfo$gene_id %in% savegenenames){
      emis1_complete  <-
        (as.numeric(Fp1[[as.character(singlegeneinfo[1, 1])]]))[c(singlegeneinfo[1, 2]:singlegeneinfo[1, 3])]
      emis2_complete  <-
        (as.numeric(Fp2[[as.character(singlegeneinfo[1, 1])]]))[c(singlegeneinfo[1, 2]:singlegeneinfo[1, 3])]
    }
    
    
  }else{
    start <- singlegeneinfo[1, 2] + endshorten
    end   <- singlegeneinfo[1, 3] - startshorten
    emis1  <-
      rev((as.integer(Fm1[[as.character(singlegeneinfo[1, 1])]]))
          [c(start:end)])
    #Note here the minus strand has been reversed!
    emis2  <-
      rev((as.integer(Fm2[[as.character(singlegeneinfo[1, 1])]]))
          [c(start:end)])
    #Note here the minus strand has been reversed!
    
    if(singlegeneinfo$gene_id %in% savegenenames){
      emis1_complete  <-
        rev((as.integer(Fm1[[as.character(singlegeneinfo[1, 1])]]))
            [c(singlegeneinfo[1, 2]:singlegeneinfo[1, 3])])
      #Note here the minus strand has been reversed!
      emis2_complete  <-
        rev((as.integer(Fm2[[as.character(singlegeneinfo[1, 1])]]))
            [c(singlegeneinfo[1, 2]:singlegeneinfo[1, 3])])
      #Note here the minus strand has been reversed!
    }
    
  }
  
  genedepth1 <- sum(emis1)
  genedepth2 <- sum(emis2)
  
  if(genedepth1 == 0 | genedepth2 == 0){
    return(NULL)
  }
  
  
  emis1_adj <- emis1/(depth1/10^6)
  emis2_adj <- emis2/(depth2/10^6)
  
  if(singlegeneinfo$gene_id %in% savegenenames){
    emis1_complete <- emis1_complete/(depth1/10^6)
    emis2_complete <- emis2_complete/(depth2/10^6)
    emis1_complete <- S4Vectors::Rle(emis1_complete)
    emis2_complete <- S4Vectors::Rle(emis2_complete)
  }else{
    emis1_complete <- NULL
    emis2_complete <- NULL
  }
  
  
  genedepth1_adj <- sum(emis1_adj)
  genedepth2_adj <- sum(emis2_adj)
  
  weight1 <- genedepth1_adj/mean(c(genedepth1_adj, genedepth2_adj))
  weight2 <- genedepth2_adj/mean(c(genedepth1_adj, genedepth2_adj))
  
  emis1_adj <- emis1_adj/weight1
  emis2_adj <- emis2_adj/weight2
  
  
  #Add pseudo count 1!!!
  emis1_adj <- emis1_adj + 1
  emis2_adj <- emis2_adj + 1
  
  #Divide the genebody into 40 bins and calculate the read count of each bin
  size <- floor(length(emis1_adj)/window_num)
  to <- (length(emis1_adj) %/% size)*size
  starts <- seq(1, to, size)
  genecov1 <- S4Vectors::Rle(emis1_adj)
  genecov2 <- S4Vectors::Rle(emis2_adj)
  vi1 <- XVector::Views(genecov1, start = starts, width = size)
  vi2 <- XVector::Views(genecov2, start = starts, width = size)
  windows1 <- XVector::ViewSums(vi1)
  windows2 <- XVector::ViewSums(vi2)
  
  #Generate the bin ratios between time points
  gene <- windows2/windows1
  
  if(sum(is.finite(gene)) <= length(gene)/2){
    
    return(NULL)
    
  }
  
  
  gene <- data.frame(ratio_adj = gene, stringsAsFactors = FALSE)
  
  res <- list(ratio_adj = gene, windowsize = size, emis1_adj = emis1_adj, emis2_adj = emis2_adj,
              emis1_complete = emis1_complete, emis2_complete = emis2_complete)
  
  return(res)
  
}

expandratio <- function(expand1 = binres$emis1_adj[startcoord:endcoord],
                        expand2 = binres$emis2_adj[startcoord:endcoord]){
  
  expand_ratio <- expand2/expand1
  
  if(sum(is.finite(expand_ratio)) <= length(expand_ratio)/3){
    
    return(NULL)
    
  }
  
  if(sum(is.finite(expand_ratio)) < length(expand_ratio)){
    
    expand_ratio <- (expand2 + 10^-5)/(expand1 + 10^-5)
    
  }
  
  expand_ratio <- data.frame(ratio_adj = expand_ratio, stringsAsFactors = FALSE)
  
  return(expand_ratio)
  
}

transpoint <- function(ratios){
  
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

plotpairs <- function(reads1 = as.numeric(genefpm1), reads2 = as.numeric(genefpm2),
                      endtype = c('TSS', 'TTS'), distance,
                      genesym = geneinfo$gene_id, genelen = geneinfo$width, timediff = time,
                      saveplot = FALSE){
  
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
  
  
  reads1$condition <- 'time2/time1'
  reads1$xcoord <- as.numeric(row.names(reads1))
  reads_totals <- reads1
  ylabel <- 'Ratio'
  
  if(!is.null(reads2)){
    reads1$condition <- 'time1'
    reads2$condition <- 'time2'
    reads2$xcoord <- as.numeric(row.names(reads2))
    reads_totals <- rbind(reads1, reads2)
    ylabel <- 'FPM'
    reads_totals$condition <- factor(reads_totals$condition, levels = c('time2', 'time1'),
                                     ordered = TRUE)
  }
  
  names(reads_totals)[1] <- 'reads'
  
  plottitle <- paste0('(Distance = ', distance, ', Length = ', genelen, ', Time difference = ', timediff, 'min)')
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
    ggplot2::theme(panel.spacing = grid::unit(0, 'lines'))
  
  if(saveplot == FALSE){
    print(p)
  }
  
  if(saveplot == TRUE){
    return(p)
  }
  
  
  
}

inferfunction <- function(datalist = parallellist[[1]], depth1 = reads1depth, depth2 = reads2depth,
                          savegenesyms = savegenenames, plotgenesyms = plotgenenames, time = time){
  
  
  ###################
  #Other default parameters
  window_num <- 40
  ###################
  #Convert the 2 bams file to strand separated wig files
  
  strandsign <- as.character(datalist$subgenes@strand@values)
  
  if('+' %in% strandsign){
    
    Fp1 <- towig(reads = datalist$subreads1, strand = '+')
    Fp2 <- towig(reads = datalist$subreads2, strand = '+')
    
  }
  
  if('-' %in% strandsign){
    
    Fm1 <- towig(reads = datalist$subreads1, strand = '-')
    Fm2 <- towig(reads = datalist$subreads2, strand = '-')
    
  }
  
  if(!('+' %in% strandsign)){
    
    Fp1 <- NULL
    Fp2 <- NULL
    
  }
  
  if(!('-' %in% strandsign)){
    
    Fm1 <- NULL
    Fm2 <- NULL
    
  }
  
  
  ###################
  startshorten <- 1000
  endshorten <- 1000
  
  gene_ids <- c()
  distances <- c()
  wilps1 <- c()
  wilps2 <- c()
  activemeans1 <- c()
  inhibitmeans1 <- c()
  activemeans2 <- c()
  inhibitmeans2 <- c()
  
  genes <- datalist$subgenes
  
  singlegenereads1 <- list()
  singlegenereads2 <- list()
  
  binresplotlist <- list()
  expandplotlist <- list()
  i <- 1
  for(i in 1:length(genes)){
    ###################
    #print(i)
    gene_len <- genes@ranges@width[i]
    
    if(gene_len < 2*(startshorten + endshorten)){
      next()
    }
    
    geneinfo <- as.data.frame(genes[i])
    binres <- binratio(singlegeneinfo = geneinfo,
                       Fp1 = Fp1, Fp2 = Fp2, Fm1 = Fm1, Fm2 = Fm2, window_num = window_num,
                       startshorten = startshorten, endshorten = endshorten,
                       depth1 = depth1, depth2 = depth2, savegenenames = savegenesyms)
    
    if(is.null(binres)){
      next()
    }
    
    ###################
    #Begin to infer
    
    point_idx <- transpoint(ratios = binres$ratio_adj)
    
    if(is.null(point_idx)){
      next()
    }
    
    if((geneinfo$gene_id %in% savegenesyms) & (plotgenesyms == TRUE)){
      binresplot <- plotpairs(reads1 = binres$ratio_adj$ratio_adj, reads2 = NULL,
                              endtype = c(paste0('TSS+', startshorten), paste0('TTS-', endshorten)),
                              distance = point_idx$point, genesym = geneinfo$gene_id,
                              genelen = length(binres$ratio_adj$ratio_adj), timediff = time,
                              saveplot = TRUE)
      binresplotlist[[geneinfo$gene_id]] <- binresplot
    }
    
    ###################
    forward <- round(binres$windowsize * (point_idx$point - 1))
    backward <- round(binres$windowsize * (point_idx$point + 1))
    
    #focus on the single window with transition point
    
    startcoord <- 1 + forward
    endcoord <- backward
    
    expand_ratio <- expandratio(expand1 = binres$emis1_adj[startcoord:endcoord],
                                expand2 = binres$emis2_adj[startcoord:endcoord])
    
    
    if(is.null(expand_ratio)){
      final_point <- list(point = binres$windowsize, front = NA, latter = NA, wilp = NA)
    }else{
      
      final_point <- transpoint(ratios = expand_ratio)
      
    }
    
    if((geneinfo$gene_id %in% savegenesyms) & (plotgenesyms == TRUE)){
      
      expandplot <- plotpairs(reads1 = expand_ratio$ratio_adj, reads2 = NULL,
                              endtype = c(paste0('TSS+', forward), paste0('TSS+', (backward - 1))),
                              distance = final_point$point, genesym = geneinfo$gene_id,
                              genelen = length(expand_ratio$ratio_adj), timediff = time,
                              saveplot = TRUE)
      expandplotlist[[geneinfo$gene_id]] <- expandplot
      
    }
    
    distance <- startshorten + forward + final_point$point
    
    gene_ids <- c(gene_ids, geneinfo[1, 6])
    distances <- c(distances, distance)
    
    wilps1 <- c(wilps1, point_idx$wilp)
    wilps2 <- c(wilps2, final_point$wilp)
    
    activemeans1 <- c(activemeans1, point_idx$front)
    inhibitmeans1 <- c(inhibitmeans1, point_idx$latter)
    
    activemeans2 <- c(activemeans2, final_point$front)
    inhibitmeans2 <- c(inhibitmeans2, final_point$latter)
    
    singlegenereads1[[geneinfo$gene_id]] <- binres$emis1_complete
    singlegenereads2[[geneinfo$gene_id]] <- binres$emis2_complete
    
    
  }
  
  subreport <- data.frame(gene_id = gene_ids,
                          distance = distances,
                          frontbinratio = activemeans1,
                          latterbinratio = inhibitmeans1,
                          diffbinratio = inhibitmeans1 - activemeans1,
                          binpval = wilps1,
                          frontextendratio = activemeans2,
                          latterextendratio = inhibitmeans2,
                          diffextendratio = inhibitmeans2 - activemeans2,
                          extendpval = wilps2,
                          stringsAsFactors = FALSE)
  
  reslist <- list(subreport = subreport,
                  singlegenereads1 = singlegenereads1,
                  singlegenereads2 = singlegenereads2,
                  binresplotlist = binresplotlist,
                  expandplotlist = expandplotlist)
  
  return(reslist)
  
}



#'Calculate gene elongation rate from a pair of Pro-seq or Gro-seq data
#'
#'Calculate gene elongation rate from a pair of Pro-seq data, or Gro-seq data, 
#'using an LSS (least sum of square) based method.
#'
#'@param time1file The reference Pro-seq/Gro-seq bam file name, corresponding 
#'  to 0min. Can be a string indicating the directory of the file, or a 
#'  GAlignmentPairs object, or a GRanges object from the original bam file.
#'@param time2file The treatment Pro-seq/Gro-seq bam file name, corresponding 
#'  to a specific time point (e.g. 15min if DRB was used for 15min). Can be a 
#'  string indicating the directory of the file, or a GAlignmentPairs object, 
#'  or a GRanges obejct from the original bam file.
#'@param targetfile The genes whose transcriptional elongation rate need to be 
#'  analyzed. If it is NULL, all the genes longer than at least 4000bp in the 
#'  genome specified by the parameter \code{genomename} will be analyzed. If 
#'  provied by the user, columns named as chr, start, end, strand, and gene_id 
#'  are required. All genes should have a length greater than at least 4000bp.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetfile} is NULL.
#'@param time The time difference between \code{time1file} and 
#'  \code{time2file}, using min as its unit.
#'@param strandmethod The strand specific method used when preparing the 
#'  sequencing library, can be 1 for directional ligation method and 2 for 
#'  dUTP method. If the sample is sequenced using a single strand method, set 
#'  it as 0.
#'@param threads Number of threads to do the parallelization. Default is 1.
#'@param savegenenames For which genes their concrete rate inference results 
#'  need to be saved, including the binned ratio values between their data in 
#'  \code{timefile2} and \code{timefile1}, the extended ratio values for 
#'  specific bins, etc.
#'@param plotgenenames Whether to plot the binned and extended ratio values 
#'  for the genes provided by \code{savegenenames}.
#'@param genelencutoff The cutoff on gene length (bp). Only genes longer than 
#'  this cutoff will be considerred. If its value is NULL or less than 4000, 
#'  genes longer than 4000bp will be considerred. Default is 70000.
#'@param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM 
#'  value greater than the cutoff in \code{time1file} will be considerred. 
#'  Default is 1.
#'@return A list including a data.frame recording the transcription elongation 
#'  rates inferred using the LSS method, as well as other information, such as 
#'  the significance of results.
#'@export
calrate <- function(time1file, time2file,
                    targetfile, genomename = NULL,
                    time,
                    strandmethod = 1, 
                    threads = 1,
                    savegenenames, plotgenenames = TRUE,
                    genelencutoff = 70000,
                    fpkmcutoff = 1){
  
  
  #Libraries required
  #library(GenomicAlignments)
  
  ##################
  #Prepare the genes to GRanges object
  startshorten <- 1000
  endshorten <- 1000
  
  if(!is.null(targetfile)){
    
    if(file.exists(targetfile)){
      
      genes <- read.table(targetfile, sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE, quote = '',
                          check.names = FALSE)
      
      genes$start <- genes$start + 1
      
      genes <- genes[(abs(genes$end - genes$start) + 1) >= max(2*(startshorten + endshorten), genelencutoff),]
      
      geneframe <- genes
      
      genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      #tstframe <- data.frame(seqnames = c('chr1', 'chr2'), start = c(1, 2500), end = c(100000, 225000),
      #                       strand = c('+', '-'), gene_id = c('tst1', 'tst2'), stringsAsFactors = FALSE)
      #geneframe <- rbind(geneframe[1:3,], tstframe, geneframe[4:47,])
      
      genes <- GenomicRanges::GRanges(seqnames = genes$seqnames,
                                      ranges = IRanges::IRanges(start = genes$start, end = genes$end),
                                      strand = genes$strand,
                                      gene_id = genes$gene_id)
      genes <- GenomicAlignments::sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- get(paste0(genomename, '.packagegenes'))
    genes <- GenomicRanges::GRanges(genes)
    
    genes <- genes[genes@ranges@width >= max(2*(startshorten + endshorten), genelencutoff)]
    
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- GenomeInfoDb::keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
                                         pruning.mode = 'coarse')
    
    genes$updis <- NULL
    genes$dndis <- NULL
    
    geneframe <- as.data.frame(genes)
    
    for(i in 1:ncol(geneframe)){
      
      if(is.factor(geneframe[,i])){
        geneframe[,i] <- as.character(geneframe[,i])
      }
      
    }
    
  }
  
  #genes <- subset(genes, gene_id %in% savegenenames)
  
  
  
  ###################
  #Read in the strand-specific paired-end data
  
  reads1 <- parsereadsfile(readsfile = time1file, strandmethod = strandmethod)
  reads2 <- parsereadsfile(readsfile = time2file, strandmethod = strandmethod)
  
  reads1depth <- length(reads1)
  reads2depth <- length(reads2)
  
  #FPKM screen###
  if(!is.null(fpkmcutoff)){
    
    genes <- getfpkm(features = genes, reads = reads1, readsdepth = reads1depth)
    genes <- genes[genes$fpkm > fpkmcutoff]
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
  }
  
  #Prepare for parallel######
  if(threads > length(genes)){
    threads <- length(genes)
  }
  
  subsize <- ceiling(length(genes)/threads)
  
  parallellist <- list()
  t <- 1
  for(t in 1:threads){
    subgenes <- genes[((t - 1)*subsize + 1):min(t*subsize, length(genes)),]
    
    subgeneranges <- range(subgenes)
    
    subreads1 <- IRanges::subsetByOverlaps(reads1, subgeneranges)
    subreads2 <- IRanges::subsetByOverlaps(reads2, subgeneranges)
    
    subreads <- c(subreads1, subreads2)
    subreadsranges <- range(subreads)
    subgenes <- IRanges::subsetByOverlaps(subgenes, subreadsranges)
    rm(subreads)
    
    subgenes <- GenomicAlignments::sort(subgenes)
    subreads1 <- GenomicAlignments::sort(subreads1)
    subreads2 <- GenomicAlignments::sort(subreads2)
    
    if(length(subgenes) == 0 | length(subreads1) == 0 | length(subreads2) == 0){
      next()
    }
    
    unitlist <- list(subgenes = subgenes, subreads1 = subreads1, subreads2 = subreads2)
    parallellist[[t]] <- unitlist
    
  }
  
  if(length(parallellist) == 0){
    return(NULL)
  }
  
  #Use parallel#####
  #date()
  #subreports <- lapply(X = parallellist, FUN = inferfunction, depth1 = reads1depth, depth2 = reads2depth,
  #                     savegenesyms = savegenenames, plotgenesyms = plotgenenames, time = time)
  #date()
  
  if(threads == 1){
    
    subreports <- list()
    
    subreports[[1]] <- inferfunction(datalist = parallellist[[1]],
                                     depth1 = reads1depth, depth2 = reads2depth,
                                     savegenesyms = savegenenames,
                                     plotgenesyms = plotgenenames,
                                     time = time)
    
    
  }else{
    
    #library(doParallel)
    #threads <- 8
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    
    #print(paste0('Actual core num is ', min(threads, cores)))
    
    doParallel::registerDoParallel(cl)
    
    date()
    `%dopar%` <- foreach::`%dopar%`
    subreports <- foreach::foreach(datalist = parallellist,
                                   .export = ls(name = globalenv())) %dopar% {inferfunction(datalist,
                                                                                            depth1 = reads1depth, depth2 = reads2depth,
                                                                                            savegenesyms = savegenenames,
                                                                                            plotgenesyms = plotgenenames,
                                                                                            time = time)}
    #.export character vector of variables to export. This can be useful when accessing a variable that isn't
    #defined in the current environment. The default valule is NULL
    date()
    
    parallel::stopCluster(cl)
    
  }
  
  
  for(i in 1:length(subreports)){
    if(i == 1){
      report <- subreports[[i]]$subreport
      genereads1 <- subreports[[i]]$singlegenereads1
      genereads2 <- subreports[[i]]$singlegenereads2
      binplots <- subreports[[i]]$binresplotlist
      expandplots <- subreports[[i]]$expandplotlist
    }else{
      report <- rbind(report, subreports[[i]]$subreport)
      genereads1 <- c(genereads1, subreports[[i]]$singlegenereads1)
      genereads2 <- c(genereads2, subreports[[i]]$singlegenereads2)
      binplots <- c(binplots, subreports[[i]]$binresplotlist)
      expandplots <- c(expandplots, subreports[[i]]$expandplotlist)
    }
  }
  
  report$time <- time
  report$rate <- report$distance/report$time
  
  report$binpadj <- p.adjust(p = report$binpval, method = 'BH')
  report$extendpadj <- p.adjust(p = report$extendpval, method = 'BH')
  
  report$significance <- 'nonsignificant'
  
  report$significance[report$diffbinratio > 0 & report$binpadj < 0.05] <- 'significant'
  
  
  geneframe <- geneframe[c('gene_id', 'seqnames', 'start', 'end', 'strand',
                           setdiff(names(geneframe), c('gene_id', 'seqnames', 'start', 'end', 'strand')))]
  
  names(report)[1] <- 'gene_id'
  
  
  sigreport <- subset(report, significance == 'significant')
  nonsigreport <- subset(report, significance == 'nonsignificant')
  
  sigreport <- sigreport[order(sigreport$binpadj, sigreport$binpval,
                               sigreport$extendpadj, sigreport$extendpval,
                               -abs(sigreport$diffbinratio), -abs(sigreport$diffextendratio),
                               -sigreport$rate),]
  nonsigreport <- nonsigreport[order(nonsigreport$binpadj, nonsigreport$binpval,
                                     nonsigreport$extendpadj, nonsigreport$extendpval,
                                     -abs(nonsigreport$diffextendratio), -abs(nonsigreport$diffextendratio),
                                     -nonsigreport$rate),]
  report <- rbind(sigreport, nonsigreport)
  row.names(report) <- 1:nrow(report)
  
  geneframe <- subset(geneframe, gene_id %in% report$gene_id)
  
  geneframe <- geneframe[match(report$gene_id, geneframe$gene_id),]
  row.names(geneframe) <- 1:nrow(geneframe)
  report <- cbind(geneframe, report[-1])
  
  names(report)[2] <- 'chr'
  
  maincols <- c('gene_id',
                'distance', 'time', 'rate',
                'significance', 'binpadj', 'binpval',
                'frontbinratio', 'latterbinratio', 'diffbinratio',
                'chr', 'start', 'end', 'strand',
                'extendpadj', 'extendpval',
                'frontextendratio', 'latterextendratio', 'diffextendratio')
  
  othercols <- setdiff(names(report), maincols)
  allcols <- c(maincols, othercols)
  report <- report[allcols]
  row.names(report) <- 1:nrow(report)
  report$genewidth <- report$end - report$start + 1
  
  res <- list()
  res[['report']] <- report
  res[['genereads1']] <- genereads1
  res[['genereads2']] <- genereads2
  
  
  if(!is.null(savegenenames) & plotgenenames == TRUE){
    for(i in 1:length(savegenenames)){
      savegenename <- savegenenames[i]
      
      print(binplots[[savegenename]])
      print(expandplots[[savegenename]])
      
      genefpm1 <- genereads1[[savegenename]]
      genefpm2 <- genereads2[[savegenename]]
      
      
      plotpairs(reads1 = as.numeric(genefpm1), reads2 = as.numeric(genefpm2),
                endtype = c('TSS', 'TTS'), distance = report$distance[report$gene_id == savegenename],
                genesym = savegenename, genelen = report$genewidth[report$gene_id == savegenename],
                timediff = time)
      
    }
    
  }
  
  
  return(res)
  
  
}



#'Calculate gene elongation rate for multiple pairs of Pro-seq or Gro-seq data
#'
#'Calculate gene elongation rate for multiple pairs of Pro-seq data, or 
#'Gro-seq data, using an LSS (least sum of square) based method.
#'
#'@param time1files The reference Pro-seq/Gro-seq bam file names, 
#'  corresponding to 0min. Should be a vector with elements as strings 
#'  indicating the directores of the files.
#'@param time2files The treatment Pro-seq/Gro-seq bam file names, 
#'  corresponding to several specific time points (e.g. 15min if DRB was used 
#'  for 15min). Should be a vector with elements as strings indicating the 
#'  directores of the files.
#'@param targetfile The genes whose transcriptional elongation rate need to be 
#'  analyzed. If it is NULL, all the genes longer than at least 4000bp in the 
#'  genome specified by the parameter \code{genomename} will be analyzed. If 
#'  provied by the user, columns named as chr, start, end, strand, and gene_id 
#'  are required. All genes should have a length greater than at least 4000bp.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetfile} is NULL.
#'@param time The time differences between time1files and time2files, using 
#'  min as its unit. Should be a vector with each element as the time 
#'  difference between the time1file and time2file presented in the 
#'  corresponding positions in \code{time1files} vector and \code{time2files} 
#'  vector.
#'@param strandmethod The strand specific method used when preparing the 
#'  sequencing library, can be 1 for directional ligation method and 2 for 
#'  dUTP method. If the sample is sequenced using a single strand method, set 
#'  it as 0.
#'@param threads Number of threads to do the parallelization. Default is 1.
#'@param mergerefs Whether merge the data of all time1 reference files 
#'  together to one reference, and then use it as a uniform reference for all 
#'  \code{time2files}. Default is TRUE.
#'@param mergecases Whether merge the data of all time2 case files together to 
#'  one case file. Default is FALSE.
#'@param savegenenames For which genes their concrete rate inference results 
#'  need to be saved, including the binned ratio values between their data in 
#'  \code{time2files} and \code{time1files}, the extended ratio values for 
#'  specific bins, etc.
#'@param plotgenenames Whether to plot the binned and extended ratio values 
#'  for the genes provided by \code{savegenenames}
#'@param genelencutoff The cutoff on gene length (bp). Only genes longer than 
#'  this cutoff will be considerred. If its value is NULL or less than 4000, 
#'  genes longer than 4000bp will be considerred. Default is 70000.
#'@param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM 
#'  value greater than the cutoff in time1file will be considerred. Default is 
#'  1.
#'@return A list with each of its element as a sub-list recording the 
#'  transcription elongation rate inferred for each pair of Pro-seq or Gro-seq 
#'  data using the LSS method, as well as other information, such as the 
#'  significance of result.
#'@export
mcalrate <- function(time1files,
                     time2files,
                     targetfile,
                     genomename = NULL,
                     times,
                     strandmethod = 1, 
                     threads = 1,
                     mergerefs = TRUE, mergecases = FALSE,
                     savegenenames, plotgenenames = TRUE,
                     genelencutoff = 70000,
                     fpkmcutoff = 1){
  
  time1files <- getreadslist(timefiles = time1files, mergefiles = mergerefs, strandmode = strandmethod)
  time2files <- getreadslist(timefiles = time2files, mergefiles = mergecases, strandmode = strandmethod)
  
  timefilelength <- length(time2files)
  
  reslist <- list()
  i <- 1
  for(i in 1:timefilelength){
    
    time2file <- time2files[[i]]
    time1file <- time1files[[min(length(time1files), i)]]
    time <- times[i]
    
    res <- calrate(time1file = time1file, time2file = time2file,
                   targetfile = targetfile, genomename = genomename,
                   time = time,
                   strandmethod = strandmethod, 
                   threads = threads,
                   savegenenames = savegenenames, plotgenenames = plotgenenames,
                   genelencutoff = genelencutoff,
                   fpkmcutoff = fpkmcutoff)
    
    reslist[[i]] <- res
    
  }
  
  return(reslist)
  
}


