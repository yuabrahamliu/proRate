
#Transition point inference#####
removeoverlap <- function(coords = genecoords, ignorestrand = TRUE){
  
  dis <- distanceToNearest(coords, ignore.strand = ignorestrand)
  disvec <- dis@elementMetadata$distance
  
  if(sum(disvec == 0)  == 0){
    return(coords)
  }else{
    newgenecoords <- coords[disvec != 0]
    removeoverlap(coords = newgenecoords, ignorestrand = ignorestrand)
  }
  
}

geneintervalscreen <- function(origenes = genes, ignorestrand = TRUE){
  
  library(GenomicRanges)
  
  if(ignorestrand == TRUE){
    ups <- follow(origenes, unstrand(origenes))
    dns <- precede(origenes, unstrand(origenes))
  }else{
    ups <- follow(origenes, origenes)
    dns <- precede(origenes, origenes)
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
  
  strandupdis <- distance(x = origenes[!is.na(strandups)], y = origenes[strandups[!is.na(strandups)]],
                          ignore.strand = ignorestrand)
  stranddndis <- distance(x = origenes[!is.na(stranddns)], y = origenes[stranddns[!is.na(stranddns)]],
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
  
  library(GenomicRanges)
  library(Biostrings)
  
  if(genomename == 'mm10'){
    library("BSgenome.Mmusculus.UCSC.mm10")
    seqs <- getSeq(Mmusculus, origenes)
  }else if(genomename == 'hg38'){
    library("BSgenome.Hsapiens.UCSC.hg38")
    seqs <- getSeq(Hsapiens, origenes)
  }
  
  #counts_genes <- oligonucleotideFrequency(seqs, width = 6, step = 1)
  #counts_genes <- as.data.frame(counts_genes)
  #freq_genes <- colSums(counts_genes)
  #freq_genes <- freq_genes[order(-freq_genes)]
  
  seqlist <- unlist(strsplit(toString(seqs), split = ', '))
  seqframe <- data.frame(genesym = origenes$gene_id,
                         seq = seqlist, stringsAsFactors = FALSE)
  
  library(seqinr)
  calGC <- function(line){
    seqstring <- as.vector(unlist(line[2]))
    seqvector <- strsplit(seqstring, split = '')[[1]]
    gccontent <- GC(seqvector)
    return(gccontent)
  }
  
  genegc <- apply(seqframe, 1, calGC)
  genegc <- data.frame(genesym = seqframe$genesym, seqgc = genegc, stringsAsFactors = FALSE)
  
  origenes$GC <- genegc$seqgc
  
  return(origenes)
  
}


calexonfreq <- function(genelineidx, origenes, exoncoords){
  
  geneline <- origenes[genelineidx]
  geneexons <- intersect(exoncoords, geneline)
  
  if(length(geneexons) == 0){
    return(0)
  }else{
    
    geneexons <- reduce(geneexons)
    geneexons$gene_id <- geneline$gene_id
    exonlens <- sum(geneexons@ranges@width)
    genelen <- geneline@ranges@width
    exonfreq <- exonlens/genelen
    
    return(exonfreq)
    
  }
  
}

geneexonfreq <- function(origenes = genes, genomename = 'mm10'){
  
  #origenes <- readRDS(paste0(genomename, '.packagegenes.rds'))[1:100]
  
  
  library(AnnotationDbi)
  
  if(genomename == 'hg38'){
    
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    
    exoncoords <- exons(TxDb.Hsapiens.UCSC.hg38.knownGene)
    
  }else if(genomename == 'mm10'){
    
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    
    exoncoords <- exons(TxDb.Mmusculus.UCSC.mm10.knownGene)
    
  }
  
  geneexonfreqs <- sapply(X = 1:length(origenes), FUN = calexonfreq, 
                          origenes = origenes, exoncoords = exoncoords)
  
  origenes$exon <- geneexonfreqs
  
  return(origenes)
  
}


makegenelist <- function(genomename = 'mm10', save = TRUE){
  
  "
  Generate a GRanges object of UCSC known genes in sepcific genome.
  
  Generate a GRanges object of UCSC known genes in specific genome.
  
  @param genomename The name of the genome, can be 'mm10' for mouse or 'hg38' for human
  @param save Whether save the result GRanges as a file in work directory, default is TRUE
  @return A GRanges object of UCSC known genes, with several gene information provided, 
  including seqnames, ranges, strands, gene symbols, distance to the nearest upstream 
  gene, distance to the nearest downstream gene, gene GC content, and gene exon portion. 
  Genes overlapping with each other will not be included.
  @export
  "
  
  library(AnnotationDbi)
  
  if(genomename == 'hg38'){
    
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    
    genecoords <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
    #Select genes with unique gene symbols
    gene_ids <- genecoords$gene_id
    genesyms <- AnnotationDbi::select(x = org.Hs.eg.db, keys = gene_ids,
                                      columns = 'SYMBOL', keytype = 'ENTREZID')
    
    targetchrs <- paste0('chr', c(seq(1, 22, 1), 'X', 'Y', 'M'))
    
  }else if(genomename == 'mm10'){
    
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    
    genecoords <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
    #Select genes with unique gene symbols
    gene_ids <- genecoords$gene_id
    genesyms <- AnnotationDbi::select(x = org.Mm.eg.db, keys = gene_ids,
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
  genecoords <- keepSeqlevels(x = genecoords, value = targetchrs, pruning.mode = 'coarse')
  
  #Remove overlapped genes, no matter the strands are different or not###
  genecoords <- removeoverlap(coords = genecoords, ignorestrand = TRUE)
  
  #Add the information of the upstream distance and downstream distance to the nearest genes,
  #no matter the strands are different or not###
  genecoords <- geneintervalscreen(origenes = genecoords, ignorestrand = TRUE)
  
  
  genecoords <- keepSeqlevels(x = genecoords, value = unique(as.character(genecoords@seqnames)),
                              pruning.mode = 'coarse')
  
  genecoords <- genegccont(origenes = genecoords, genomename = genomename)
  genecoords <- geneexonfreq(origenes = genecoords, genomename = genomename)
  
  genecoords <- sort(genecoords)
  genes <- genecoords
  names(genes) <- 1:length(genes)
  
  if(save == TRUE){
    saveRDS(genes, paste0(genomename, '.packagegenes.rds'))
  }
  
  return(genes)
  
}

getfpkm <- function(features = genes, reads = reads1,
                    readsdepth = reads1depth){
  
  fcounts <- summarizeOverlaps(features = features, reads = reads,
                               mode = 'IntersectionNotEmpty',
                               singleEnd = FALSE, ignore.strand = FALSE)
  #Add singleEnd = FALSE for pair-end
  #Add ignore.strand = FALSE for strand specific
  
  fcounts <- as.data.frame(assays(fcounts)$counts)
  featurewidth <- features@ranges@width
  
  factorM <- readsdepth/10^6
  factorK <- featurewidth/10^3
  
  fpkms <- fcounts/factorK/factorM
  features$fpkm <- fpkms$reads
  
  return(features)
  
}

parsereadsfile <- function(readsfile = time1file, strandmethod = 1){
  
  library(GenomicAlignments)
  
  if(is.character(readsfile)){
    if(file.exists(readsfile)){
      
      if(strandmethod != 0){
        reads_ori <- readGAlignmentPairs(readsfile, strandMode = strandmethod)
      }else{
        reads_ori <- readGAlignments(readsfile)
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
    reads <- reads[seqnames(reads) == chrom,]
  }
  
  
  readsList <- split(reads, seqnames(reads))
  #important!
  #'split' divides the data in a vector-like object 'x' into the
  #groups defined by 'f'.
  #Here split the bam file data into different chroms
  
  # Change reads strand
  readsList <- endoapply(readsList, function(x){
    
    if(strand == "*"){
      strand(x) <- "*"
    }else{
      x <- x[strand(x) == strand,]
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
    seqlevels(x) <- seqlevelsInUse(x)
    #seqlevels contains all the chrs in a factor
    #seqlevelsInUse only contains the chrs present in x
    
    cov <- coverage(x)[[1]]
    #Important
    #Caluculate the coverage per base, like wig file, per chrom
    to <- length(cov)
    #length(cov) = seqlengths(x), it is chrom length
    starts <- 1:to
    vi <- Views(cov, start=starts, width=1)
    
    H[[i]] <- Rle(viewSums(vi))
    #'viewSums' calculate sums of the views in a Views or ViewsList object.
    #Here the total read num of each window can be calculated
    #'Rle(values, lengths)': This constructor creates an Rle instance out of
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
    emis1_complete <- Rle(emis1_complete)
    emis2_complete <- Rle(emis2_complete)
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
  genecov1 <- Rle(emis1_adj)
  genecov2 <- Rle(emis2_adj)
  vi1 <- Views(genecov1, start = starts, width = size)
  vi2 <- Views(genecov2, start = starts, width = size)
  windows1 <- viewSums(vi1)
  windows2 <- viewSums(vi2)
  
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
  
  
  
  
  library(ggplot2)
  library(scales)
  library(grid)
  
  p <- ggplot(reads_totals, mapping = aes(x = xcoord, y = reads))
  
  if(!is.null(reads2)){
    p <- p + geom_line(size = 1, color = hue_pal()(1))
  }else{
    p <- p + geom_point(size = 2, color = hue_pal()(1))
  }
  
  p <- p + scale_x_continuous(breaks = points,
                              labels = pointlabels) +
    xlab('') + ylab(ylabel) +
    geom_vline(xintercept = points[2], linetype = 2, color = 'red', size = 1) +
    facet_grid(condition~., scales = 'free_y') +
    ggtitle(plottitle) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.x=element_text(hjust=1, angle = 90)) +
    theme(panel.spacing = unit(0, 'lines'))
  
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

calrate <- function(time1file = time1file, time2file = time2file,
                    targetfile = targetfile, genomename = NULL,
                    time = 15,
                    strandmethod = 1, 
                    threads = 2,
                    savegenenames = c('Tbl1x', 'Mid1', 'Lrch2'), plotgenenames = TRUE,
                    genelencutoff = 70000,
                    fpkmcutoff = 1){
  
  "
  Calculate gene elongation rate from a pair of Pro-seq or Gro-seq data.
  
  Calculate gene elongation rate from a pair of Pro-seq data, or Gro-seq data, 
  using an LSS (least sum of square) based method.
  
  @param time1file The reference Pro-seq/Gro-seq bam file name, corresponding to 0min.
                   Can be a string indicating the directory of the file, or a GAlignmentPairs 
                   object, or a GRanges object from the original bam file
  @param time2file The treatment Pro-seq/Gro-seq bam file name, corresponding to a specific 
                   time point (e.g. 15min if DRB was used for 15min). Can be a string indicating 
                   the directory of the file, or a GAlignmentPairs object, or a GRanges obejct 
                   from the original bam file
  @param targetfile The genes whose transcriptional elongation rate need to be analyzed. If it 
                    is NULL, all the genes longer than at least 4000bp in the genome specified by the 
  parameter genomename will be analyzed. If provied by the user, columns named 
                       as chr, start, end, strand, and gene_id are required. All genes should have 
                       a length greater than at least 4000bp.
  @param genomename Specify the genome of the genes to be analyzed, when the parameter targetfile 
                    is NULL.
  @param time The time difference between time1file and time2file, using min as its unit.
  @param strandmethod The strand specific method used when preparing the sequencing library,
                      can be 1 for directional ligation method and 2 for dUTP method.
                      If the sample is sequenced using a single strand method, set it as 0.
  @param threads Number of threads to do the parallelization. Default is 1.
  @param savegenenames For which genes their concrete rate inference results need to be saved, 
                       including the binned ratio values between their data in timefile2 and 
                       timefile1, the extended ratio values for specific bins, etc.
  @param plotgenenames Whether to plot the binned and extended ratio values for the genes 
                       provided by savegenenames
  @param genelencutoff The cutoff on gene length (bp). Only genes longer than this cutoff will be 
                       considerred. If its value is NULL or less than 4000, genes longer than 
                       4000bp will be considerred. Default is 70000.
  @param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM value greater than the 
                    cutoff will be considerred. Default is 1.
  
  @return A list including a DataFrame recording the transcription elongation rates inferred using 
          the LSS method, as well as other information, such as the significance of results.
  @export
  "
  
  #Libraries required
  library(GenomicAlignments)
  
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
      
      genes <- GRanges(seqnames = genes$seqnames,
                       ranges = IRanges(start = genes$start, end = genes$end),
                       strand = genes$strand,
                       gene_id = genes$gene_id)
      genes <- sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    genes <- genes[genes@ranges@width >= max(2*(startshorten + endshorten), genelencutoff)]
    
    genes <- sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
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
    genes <- sort(genes)
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
    
    subreads1 <- subsetByOverlaps(reads1, subgeneranges)
    subreads2 <- subsetByOverlaps(reads2, subgeneranges)
    
    subreads <- c(subreads1, subreads2)
    subreadsranges <- range(subreads)
    subgenes <- subsetByOverlaps(subgenes, subreadsranges)
    rm(subreads)
    
    subgenes <- sort(subgenes)
    subreads1 <- sort(subreads1)
    subreads2 <- sort(subreads2)
    
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
  
  
  library(doParallel)
  #threads <- 8
  
  cores <- detectCores()
  cl <- makeCluster(min(threads, cores))
  
  #print(paste0('Actual core num is ', min(threads, cores)))
  
  registerDoParallel(cl)
  
  date()
  subreports <- foreach(datalist = parallellist,
                        .export = ls(name = globalenv())) %dopar% {inferfunction(datalist,
                                                                                 depth1 = reads1depth, depth2 = reads2depth,
                                                                                 savegenesyms = savegenenames,
                                                                                 plotgenesyms = plotgenenames,
                                                                                 time = time)}
  #.export character vector of variables to export. This can be useful when accessing a variable that isn't
  #defined in the current environment. The default valule is NULL
  date()
  
  stopCluster(cl)
  
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

mcalrate <- function(time1files = c(wt0file, ko0file),
                     time2files = c(wt15file, wt30file),
                     targetfile = targetfile,
                     genomename = NULL,
                     times = c(15, 30),
                     strandmethod = 1, 
                     threads = 1,
                     mergerefs = TRUE, mergecases = FALSE,
                     savegenenames = c('Tbl1x', 'Mid1', 'Lrch2'), plotgenenames = TRUE,
                     genelencutoff = 70000,
                     fpkmcutoff = 1){
  
  "
  Calculate gene elongation rate from multiple pairs of Pro-seq or Gro-seq data.
  
  Calculate gene elongation rate from multiple pairs of Pro-seq data, or Gro-seq data, 
  using an LSS (least sum of square) based method.
  
  @param time1files The reference Pro-seq/Gro-seq bam file names, corresponding to 0min.
                    Should be a vector with elements as strings indicating the directores of 
                    the files.
  @param time2files The treatment Pro-seq/Gro-seq bam file names, corresponding to several specific 
                    time points (e.g. 15min if DRB was used for 15min). Should be a vector with 
                    elements as strings indicating the directores of the files.
  @param targetfile The genes whose transcriptional elongation rate need to be analyzed. 
                    If it is NULL, all the genes longer than at least 4000bp in the genome 
                    specified by the parameter genomename will be analyzed. If provied by 
                    the user, columns named as chr, start, end, strand, and gene_id are 
                    required. All genes should have a length greater than at least 4000bp.
  @param genomename Specify the genome of the genes to be analyzed, when the parameter targetfile 
                    is NULL.
  @param time The time differences between time1files and time2files, using min as its unit. 
              Should be a vector with each element as the time difference between the time1file 
              and time2file presented in the corresponding positions in time1files vector and 
              time2files vector.
  @param strandmethod The strand specific method used when preparing the sequencing library,
                      an be 1 for directional ligation method and 2 for dUTP method.
                      If the sample is sequenced using a single strand method, set it as 0.
  @param threads Number of threads to do the parallelization. Default is 1.
  @param mergerefs Whether merge the data of all time1 reference files together as one reference, 
                   and then use it as reference for all time2files. Default is TRUE.
  @param mergecases Whether merge the data of all time2 case files together as one case file. 
                    Default is FALSE.
  @param savegenenames For which genes their concrete rate inference results need to be saved, 
                       including the binned ratio values between their data in timefile2 and 
                       timefile1, the extended ratio values for specific bins, etc.
  @param plotgenenames Whether to plot the binned and extended ratio values for the genes 
                       provided by savegenenames
  @param genelencutoff The cutoff on gene length (bp). Only genes longer than this cutoff will be 
                       considerred. If its value is NULL or less than 4000, genes longer than 
                       4000bp will be considerred. Default is 70000.
  @param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM value greater than the 
                    cutoff will be considerred. Default is 1.
  
  @return A list and each of its element is a sub-list recording the transcription elongation rate 
          inferred from each pair of Pro-seq or Gro-seq data using the LSS method, as well as other 
          information, such as the significance of result.
  @export
  "
  
  time1files <- getreadslist(timefiles = time1files, mergefiles = mergerefs, strandmode = starndmethod)
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

#Statistics (Intra group comparison and inter group comparison)#####
quantileattr <- function(dat = rates, quantiles = quantilenum){
  
  absquantiles <- quantile(x = dat$absrate, probs = seq(0, 1, 1/quantiles))
  relquantiles <- quantile(x = dat$relrate, probs = seq(0, 1, 1/quantiles))
  
  absquantiles <- as.vector(absquantiles)
  relquantiles <- as.vector(relquantiles)
  
  absdat <- dat[c('gene_id', 'absrate')]
  reldat <- dat[c('gene_id', 'relrate')]
  absdat$absq <- 'quantile1'
  reldat$relq <- 'quantile1'
  
  i <- 2
  for(i in 2:length(absquantiles)){
    quantilename <- paste0('quantile', (i - 1))
    absdat$absq[absdat$absrate > absquantiles[i - 1] & absdat$absrate <= absquantiles[i]] <- quantilename
    reldat$relq[reldat$relrate > relquantiles[i - 1] & reldat$relrate <= relquantiles[i]] <- quantilename
    
  }
  
  absdat <- absdat[order(-absdat$absrate),]
  reldat <- reldat[order(-reldat$relrate),]
  row.names(absdat) <- 1:nrow(absdat)
  row.names(reldat) <- 1:nrow(reldat)
  
  resdat <- list(absres = absdat, relres = reldat)
  
  
}

fishermean <- function(line = line){
  
  a11 <- line[1]
  a12 <- line[2]
  a21 <- a22 <- mean(line)
  mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
  mat <- round(mat)
  fisherres <- fisher.test(mat)
  fisherp  <- fisherres$p.value
  
  return(fisherp)
  
}

wilcoxcomp <- function(idx, mat1 = dat1, mat2 = dat2){
  
  wilcoxres <- wilcox.test(unlist(mat1[idx,]), unlist(mat2[idx,]))
  wilcoxp <- wilcoxres$p.value
  
  return(wilcoxp)
  
}

limmacomp <- function(mat1 = dat1, mat2 = dat2){
  
  mat <- cbind(mat1, mat2)
  
  #mat <- log2(mat)
  
  pd <- data.frame(group = colnames(mat), stringsAsFactors = FALSE)
  row.names(pd) <- colnames(mat)
  pd$group <- c(rep('group1', ncol(mat1)), rep('group2', ncol(mat2)))
  
  library(limma)
  designstr <- '~ group'
  designstr <- as.formula(designstr)
  
  design <- model.matrix(designstr, data = pd)
  fit1 <- lmFit(mat, design)
  fit2 <- eBayes(fit1)
  
  allg.limma <- topTable(fit2, coef = 2, n = dim(fit1)[1])
  
  return(allg.limma)
  
}




intracompare <- function(inferres = allrateres, targetgenes = NULL, quantilenum = 4){
  
  
  if(!('report' %in% names(inferres))){
    
    i <- 1
    for(i in 1:length(inferres)){
      
      subres <- inferres[[i]]$report
      
      if(!is.null(targetgenes)){
        
        if(targetgenes == 'significant'){
          subres <- subset(subres, significance == 'significant')
        }else{
          
          subres <- subset(subres, gene_id %in% targetgenes)
          
        }
        
      }
      
      subgenes <- subres$gene_id
      
      if(i == 1){
        finalgenes <- subgenes
      }else{
        finalgenes <- intersect(finalgenes, subgenes)
      }
      
    }
    
    j <- 1
    for(j in 1:length(inferres)){
      
      subres <- inferres[[j]]$report
      subres <- subres[match(finalgenes, subres$gene_id),]
      subres$relrate <- subres$rate/(abs(subres$end - subres$start) + 1)
      subres <- subres[c('gene_id', 'rate', 'relrate')]
      names(subres) <- c('gene_id', paste0('absrate_', j), paste0('relrate_', j))
      
      subabs <- subres[,c(1, 2)]
      subrel <- subres[,c(1, 3)]
      
      if(j == 1){
        absrates <- subabs
        relrates <- subrel
      }else{
        absrates <- cbind(absrates, subabs[2])
        relrates <- cbind(relrates, subrel[2])
      }
      
      
    }
    
    absrates <- data.frame(gene_id = absrates$gene_id, absrate = rowMeans(absrates[-1]),
                           stringsAsFactors = FALSE)
    relrates <- data.frame(gene_id = relrates$gene_id, relrate = rowMeans(relrates[-1]),
                           stringsAsFactors = FALSE)
    
    rates <- cbind(absrates, relrates[-1])
    row.names(rates) <- 1:nrow(rates)
    
    
  }else{
    
    subres <- inferres$report
    
    if(!is.null(targetgenes)){
      
      if(targetgenes == 'significant'){
        subres <- subset(subres, significance == 'significant')
      }else{
        
        subres <- subset(subres, gene_id %in% targetgenes)
        
      }
      
    }
    
    subres$relrate <- subres$rate/(abs(subres$end - subres$start) + 1)
    subres <- subres[c('gene_id', 'rate', 'relrate')]
    names(subres) <- c('gene_id', 'absrate', 'relrate')
    
    rates <- subres
    row.names(rates) <- 1:nrow(rates)
    
  }
  
  res <- quantileattr(dat = rates, quantiles = quantilenum)
  
  names(res) <- c('absolute_rate', 'relative_rate')
  
  return(res)
  
}

orgintercompareinput <- function(reslist = allrateres, targetgenes = NULL){
  
  i <- 1
  for(i in 1:length(reslist)){
    
    res <- reslist[[i]]
    res <- res$report
    
    if(!is.null(targetgenes)){
      
      if(targetgenes == 'significant'){
        res <- subset(res, significance == 'significant')
      }else{
        
        res <- subset(res, gene_id %in% targetgenes)
        
      }
    }
    
    res <- res[c('gene_id', 'rate')]
    gene_id <- res$gene_id
    
    if(i == 1){
      gene_ids <- gene_id
      res <- res[match(gene_ids, res$gene_id),]
      reses <- res
    }else{
      gene_ids <- intersect(gene_ids, gene_id)
      res <- res[match(gene_ids, res$gene_id),]
      reses <- reses[match(gene_ids, reses$gene_id),]
      reses <- cbind(reses, res[-1])
    }
    
    
  }
  
  names(reses) <- c('gene_id', paste0('rate_', seq(1, (ncol(reses) - 1), 1)))
  
  return(reses)
  
}

intercompare <- function(inferresmat = rates, groupnames = c('15min', '30min')){
  
  groups <- unique(groupnames)
  if(length(groups) != 2){
    
    print('Currently only comparison between 2 groups is supported, so the data should come from 2 groups')
    
    return(NULL)
    
  }
  
  groupidces <- paste0(groupnames, '.', seq(1, length(groupnames), 1))
  names(inferresmat) <- c('gene_id', groupidces)
  
  dat1 <- groupidces[groupnames == groups[1]]
  dat2 <- groupidces[groupnames == groups[2]]
  
  dat1 <- inferresmat[c('gene_id', dat1)]
  dat2 <- inferresmat[c('gene_id', dat2)]
  
  colnames(dat1) <- c('gene_id', paste0(groups[1], '.', 1:(ncol(dat1) - 1)))
  colnames(dat2) <- c('gene_id', paste0(groups[2], '.', 1:(ncol(dat2) - 1)))
  
  
  row.names(dat1) <- dat1$gene_id
  row.names(dat2) <- dat2$gene_id
  dat1 <- dat1[2:ncol(dat1)]
  dat2 <- dat2[2:ncol(dat2)]
  
  dat1 <- log2(dat1)
  dat2 <- log2(dat2)
  
  if(ncol(dat1) < 3 & ncol(dat2) < 3){
    
    dat1 <- rowMeans(dat1)
    dat2 <- rowMeans(dat2)
    
    dat <- cbind(dat1, dat2)
    
    pvals <- apply(X = dat, MARGIN = 1, FUN = fishermean)
    padjs <- p.adjust(pvals, method = 'BH')
    log2fcs <- dat2 - dat1
    
    group1rates <- dat1
    group2rates <- dat2
    gene_ids <- names(dat1)
    
  }else{
    
    limmares <- limmacomp(mat1 = dat1, mat2 = dat2)
    limmapvals <- limmares[c('P.Value', 'adj.P.Val')]
    names(limmapvals) <- c('pval(limma)', 'padj(limma)')
    limmapvals$gene_id <- row.names(limmares)
    
    testlist <- as.list(seq(1, nrow(dat1), 1))
    
    pvals <- lapply(X = testlist, FUN = wilcoxcomp, mat1 = dat1, mat2 = dat2)
    names(pvals) <- row.names(dat1)
    pvals <- unlist(pvals)
    padjs <- p.adjust(pvals, method = 'BH')
    log2fcs <- rowMeans(dat2) - rowMeans(dat1)
    
    group1rates <- rowMeans(dat1)
    group2rates <- rowMeans(dat2)
    gene_ids <- row.names(dat1)
    
  }
  
  res <- data.frame(gene_id = gene_ids, group1rate = group1rates, group2rate = group2rates,
                    pval = pvals, padj = padjs, log2FC = log2fcs, stringsAsFactors = FALSE)
  
  names(res) <- c('gene_id', paste0(groups[1], '.rate'), paste0(groups[2], '.rate'),
                  'pval', 'padj', paste0('log2FC(', groups[2], '/', groups[1], ')'))
  row.names(res) <- 1:nrow(res)
  
  if(exists('limmapvals')){
    res <- merge(res, limmapvals, by = c('gene_id'))
  }
  
  res <- res[order(res$padj, res$pval, -abs(res[,6])),]
  row.names(res) <- 1:nrow(res)
  
  return(res)
  
}


#ssMetagene (ssMetagene)#####

revstrand <- function(oristrand){
  
  oristrand <- factor(oristrand, levels = c('+', '-'), labels = c('-', '+'))
  oristrand <- as.character(oristrand)
  return(oristrand)
  
}

convertartifical <- function(reads1 = subreads1, reads2 = subreads2, strandmethod = 1){
  
  #Change artifitial reads to true reads#####
  firststrands <- as.character(reads1@strand)
  secondstrands <- as.character(reads2@strand)
  
  if(strandmethod == 1){
    secondstrands <- revstrand(secondstrands)
  }else{
    firststrands <- revstrand(firststrands)
  }
  
  
  strand(reads1) <- firststrands
  strand(reads2) <- secondstrands
  #Use strand(reads1) to change the strand symbols, don't use reads1@strand, which will not work,
  #or rewrite reads1 <- GRanges(seqnames = reads1@seqnames, ranges = reads1@ranges,
  #strand = firststrands), which will lost other important information in the original reads1, such as
  #reads1@seqinfo, which includes chromosome length information
  
  reads <- c(reads1, reads2)
  reads <- sort(reads)
  
  return(reads)
  
}

countregion <- function(regioncount = tssfwdcount, start = 2001, width = maxtssradius + 1,
                        wholedepth = depth){
  
  end <- start + width - 1
  subcount <- regioncount[start:end]
  wholefactor <- wholedepth/10^6
  subfpkm <- sum(subcount)/wholefactor/(width/10^3)
  
  return(subfpkm)
  
}

compressgblength <- function(transfac = transfac, gene_len = gene_len, genebodylen = genebodylen,
                             counts = gbfwdcount){
  
  transfac <- floor(transfac)
  remainings <- gene_len - transfac*genebodylen
  
  part1 <- 1:(floor((gene_len - remainings*(transfac + 1))/transfac/2)*transfac)
  part2 <- (length(part1) + 1):(length(part1) + remainings*(transfac + 1))
  part3 <- (length(part1) + length(part2) + 1):gene_len
  
  part1mat <- matrix(counts[part1], ncol = transfac, byrow = TRUE)
  part2mat <- matrix(counts[part2], ncol = transfac + 1, byrow = TRUE)
  part3mat <- matrix(counts[part3], ncol = transfac, byrow = TRUE)
  
  counttrans <- c(rowMeans(part1mat), rowMeans(part2mat), rowMeans(part3mat))
  
  return(counttrans)
  
}

extendgblength <- function(transfac = transfac, gene_len = gene_len, genebodylen = genebodylen,
                           counts = gbfwdcount){
  
  revfac <- 1/transfac
  revfac <- ceiling(revfac)
  remainings <- revfac*gene_len - genebodylen
  
  part1 <- 1:(floor(remainings/2/(revfac - 1))*(revfac - 1))
  part3 <- (gene_len - (remainings - length(part1)) + 1):gene_len
  part2 <- (part1[length(part1)] + 1):(part3[1] - 1)
  
  part1trans <- rep(counts[part1], each = revfac - 1)
  part2trans <- rep(counts[part2], each = revfac)
  part3trans <- rep(counts[part3], each = revfac - 1)
  
  counttrans <- c(part1trans, part2trans, part3trans)
  
  return(counttrans)
  
  
}




mapfunction <- function(datalist = parallellist[[3]], depth = readsdepth,
                        tssradius = tssradius, ttsradius = ttsradius, genebodylen = genebodylen,
                        saveindividual = plotgenenames){
  
  #It is important to keep the seqinfo slot of datalist$subreads, because the function towig needs it to
  #align the reads in datalist$subreads to the absolute coordinates of sepcific chromosomes. If this slot
  #were not contained in datalist$subreads, the output will has wrong absolute coordinates for the reads.
  Fp <- towig(reads = datalist$subreads, strand = '+')
  Fm <- towig(reads = datalist$subreads, strand = '-')
  
  maxtssradius <- max(tssradius)
  maxttsradius <- max(ttsradius)
  genes <- datalist$subgenes
  
  chrlimits <- datalist$subreads@seqinfo@seqlengths
  names(chrlimits) <- datalist$subreads@seqinfo@seqnames
  chrlimits <- chrlimits[as.character(genes@seqnames)]
  chrlimits <- as.vector(chrlimits)
  
  starts <- genes@ranges@start
  ends <- starts + genes@ranges@width - 1
  
  tsses <- starts
  tsses[as.character(genes@strand) == '-'] <- ends[as.character(genes@strand) == '-']
  ttses <- ends
  ttses[as.character(genes@strand) == '-'] <- starts[as.character(genes@strand) == '-']
  
  tssups <- tsses - maxtssradius
  tssdns <- tsses + maxtssradius
  ttsups <- ttses - maxttsradius
  ttsdns <- ttses + maxttsradius
  
  tssups[tssups < 1] <- 1
  ttsups[ttsups < 1] <- 1
  ttsdns[ttsdns > chrlimits] <- chrlimits[ttsdns > chrlimits]
  tssdns[tssdns > chrlimits] <- chrlimits[tssdns > chrlimits]
  
  genecoords <- data.frame(seqnames = as.character(genes@seqnames), tssups = tssups, tssdns = tssdns,
                           ttsups = ttsups, ttsdns = ttsdns, strand = as.character(genes@strand),
                           starts = starts, ends = ends,
                           gene_id = as.character(genes$gene_id),
                           stringsAsFactors = FALSE)
  
  wholefactor <- depth/10^6
  
  tssfpkms <- list()
  ttsfpkms <- list()
  genebodyfpkms <- list()
  
  
  tssfwdcum <- rep(0, sum(maxtssradius, 1, maxtssradius))
  tssrevcum <- rep(0, length(tssfwdcum))
  ttsfwdcum <- rep(0, sum(maxttsradius, 1, maxttsradius))
  ttsrevcum <- rep(0, length(ttsfwdcum))
  gbfwdcum <- rep(0, genebodylen)
  gbrevcum <- rep(0, genebodylen)
  
  tssfwdindividual <- list()
  tssrevindividual <- list()
  ttsfwdindividual <- list()
  gbfwdindividual <- list()
  
  ttsrevindividual <- list()
  gbrevindividual <- list()
  
  i <- 1
  for(i in 1:nrow(genecoords)){
    #print(i)
    gene_len <- genes@ranges@width[i]
    
    if(gene_len < max(c(tssradius, ttsradius))){
      next()
    }
    
    genecoord <- genecoords[i,]
    
    genechr <- genecoord$seqnames
    genestrand <- genecoord$strand
    
    if(maxtssradius > 0){
      
      if(genestrand == '+'){
        
        tssfwdcount <- (as.numeric(Fp[[genechr]]))[(genecoord$tssups):(genecoord$tssdns)]
        tssrevcount <- (as.numeric(Fm[[genechr]]))[(genecoord$tssups):(genecoord$tssdns)]
        
      }else{
        
        tssfwdcount <- (as.numeric(Fm[[genechr]]))[(genecoord$tssups):(genecoord$tssdns)]
        tssfwdcount <- rev(tssfwdcount)
        tssrevcount <- (as.numeric(Fp[[genechr]]))[(genecoord$tssups):(genecoord$tssdns)]
        tssrevcount <- rev(tssrevcount)
        
      }
      
      j <- 1
      for(j in 1:length(tssradius)){
        
        subradius <- tssradius[j]
        subfwdfpkm <- countregion(regioncount = tssfwdcount, start = maxtssradius + 1,
                                  width = subradius + 1, wholedepth = depth)
        subrevfpkm <- countregion(regioncount = tssrevcount, start = maxtssradius - subradius + 1,
                                  width = subradius, wholedepth = depth)
        
        subfpkm <- data.frame(gene_id = genecoord$gene_id, fwdfpkm = subfwdfpkm, revfpkm = subrevfpkm,
                              stringsAsFactors = FALSE)
        names(subfpkm) <- c('gene_id', paste0('TSS', subradius, '_fwdfpkm'),
                            paste0('TSS', subradius, '_revfpkm'))
        if(j == 1){
          tsssubfpkms <- subfpkm
        }else{
          tsssubfpkms <- cbind(tsssubfpkms, subfpkm[-1])
        }
        
      }
      
      tssfpkms[[genecoord$gene_id]] <- tsssubfpkms
      
      tssfwdcum <- tssfwdcum + tssfwdcount
      tssrevcum <- tssrevcum + tssrevcount
      
      if(genecoord$gene_id %in% saveindividual){
        
        tssfwdfpm <- tssfwdcount/wholefactor
        tssrevfpm <- tssrevcount/wholefactor
        
        tssfwdindividual[[genecoord$gene_id]] <- Rle(tssfwdfpm)
        tssrevindividual[[genecoord$gene_id]] <- Rle(tssrevfpm)
        
      }
      
      
    }
    
    if(maxttsradius > 0){
      
      if(genestrand == '+'){
        
        ttsfwdcount <- (as.numeric(Fp[[genechr]]))[(genecoord$ttsups):(genecoord$ttsdns)]
        ttsrevcount <- (as.numeric(Fm[[genechr]]))[(genecoord$ttsups):(genecoord$ttsdns)]
        
      }else{
        
        ttsfwdcount <- (as.numeric(Fm[[genechr]]))[(genecoord$ttsups):(genecoord$ttsdns)]
        ttsfwdcount <- rev(ttsfwdcount)
        ttsrevcount <- (as.numeric(Fp[[genechr]]))[(genecoord$ttsups):(genecoord$ttsdns)]
        ttsrevcount <- rev(ttsrevcount)
        
      }
      
      
      j <- 1
      for(j in 1:length(ttsradius)){
        
        subradius <- ttsradius[j]
        subfwdfpkmdn <- countregion(regioncount = ttsfwdcount, start = maxttsradius + 1,
                                    width = subradius + 1, wholedepth = depth)
        subfwdfpkmup <- countregion(regioncount = ttsfwdcount, start = maxttsradius - subradius + 1,
                                    width = subradius, wholedepth = depth)
        
        subfpkm <- data.frame(gene_id = genecoord$gene_id, fwdfpkm = subfwdfpkmdn, revfpkm = subfwdfpkmup,
                              stringsAsFactors = FALSE)
        names(subfpkm) <- c('gene_id', paste0('TTS', subradius, '_fwdfpkmdn'),
                            paste0('TTS', subradius, '_fwdfpkmup'))
        if(j == 1){
          ttssubfpkms <- subfpkm
        }else{
          ttssubfpkms <- cbind(ttssubfpkms, subfpkm[-1])
        }
        
      }
      
      ttsfpkms[[genecoord$gene_id]] <- ttssubfpkms
      
      ttsfwdcum <- ttsfwdcum + ttsfwdcount
      ttsrevcum <- ttsrevcum + ttsrevcount
      
      if(genecoord$gene_id %in% saveindividual){
        
        ttsfwdfpm <- ttsfwdcount/wholefactor
        ttsrevfpm <- ttsrevcount/wholefactor
        
        ttsfwdindividual[[genecoord$gene_id]] <- Rle(ttsfwdfpm)
        ttsrevindividual[[genecoord$gene_id]] <- Rle(ttsrevfpm)
        
      }
      
      
    }
    
    if(genebodylen > 0){
      
      if(genestrand == '+'){
        
        gbfwdcount <- (as.numeric(Fp[[genechr]]))[(genecoord$starts):(genecoord$ends)]
        gbrevcount <- (as.numeric(Fm[[genechr]]))[(genecoord$starts):(genecoord$ends)]
        
      }else{
        
        gbfwdcount <- (as.numeric(Fm[[genechr]]))[(genecoord$starts):(genecoord$ends)]
        gbfwdcount <- rev(gbfwdcount)
        gbrevcount <- (as.numeric(Fp[[genechr]]))[(genecoord$starts):(genecoord$ends)]
        gbrevcount<- rev(gbrevcount)
        
      }
      
      
      gbfwdfpkm <- countregion(regioncount = gbfwdcount, start = 1, width = gene_len, wholedepth = depth)
      
      subfpkm <- data.frame(gene_id = genecoord$gene_id, fwdfpkm = gbfwdfpkm,
                            stringsAsFactors = FALSE)
      names(subfpkm) <- c('gene_id', paste0('genebody_fwdfpkmdn'))
      gbfpkms <- subfpkm
      
      genebodyfpkms[[genecoord$gene_id]] <- gbfpkms
      
      
      transfac <- gene_len/genebodylen
      
      if(transfac >= 1){
        
        gbfwdtrans <- compressgblength(transfac = transfac, gene_len = gene_len, genebodylen = genebodylen,
                                       counts = gbfwdcount)
        gbrevtrans <- compressgblength(transfac = transfac, gene_len = gene_len, genebodylen = genebodylen,
                                       counts = gbrevcount)
        
        
      }else{
        
        gbfwdtrans <- extendgblength(transfac = transfac, gene_len = gene_len, genebodylen = genebodylen,
                                     counts = gbfwdcount)
        gbrevtrans <- extendgblength(transfac = transfac, gene_len = gene_len, genebodylen = genebodylen,
                                     counts = gbrevcount)
        
      }
      
      gbfwdcum <- gbfwdcum + gbfwdtrans
      gbrevcum <- gbrevcum + gbrevtrans
      
      if(genecoord$gene_id %in% saveindividual){
        
        gbfwdfpm <- gbfwdcount/wholefactor
        gbrevfpm <- gbrevcount/wholefactor
        
        gbfwdindividual[[genecoord$gene_id]] <- Rle(gbfwdfpm)
        gbrevindividual[[genecoord$gene_id]] <- Rle(gbrevfpm)
      }
      
      
      
      
    }
    
    
  }
  
  reslist <- list()
  
  if(maxtssradius > 0){
    
    tssfwdcumfpm <- tssfwdcum/wholefactor
    tssrevcumfpm <- tssrevcum/wholefactor
    
    tssfpkms <- do.call(rbind, tssfpkms)
    row.names(tssfpkms) <- 1:nrow(tssfpkms)
    reslist[['TSSFPKMstats']] <- tssfpkms
    
    reslist[['TSSfwdFPMcum']] <- tssfwdcumfpm
    reslist[['TSSrevFPMcum']] <- tssrevcumfpm
    
    reslist[['TSSsinglefwdFPM']] <- tssfwdindividual
    reslist[['TSSsinglerevFPM']] <- tssrevindividual
    
  }
  
  if(maxttsradius > 0){
    
    ttsfwdcumfpm <- ttsfwdcum/wholefactor
    ttsrevcumfpm <- ttsrevcum/wholefactor
    
    ttsfpkms <- do.call(rbind, ttsfpkms)
    row.names(ttsfpkms) <- 1:nrow(ttsfpkms)
    reslist[['TTSFPKMstats']] <- ttsfpkms
    
    reslist[['TTSfwdFPMcum']] <- ttsfwdcumfpm
    reslist[['TTSrevFPMcum']] <- ttsrevcumfpm
    
    reslist[['TTSsinglefwdFPM']] <- ttsfwdindividual
    reslist[['TTSsinglerevFPM']] <- ttsrevindividual
    
  }
  
  if(genebodylen > 0){
    
    gbfwdcumfpm <- gbfwdcum/wholefactor
    gbrevcumfpm <- gbrevcum/wholefactor
    
    genebodyfpkms <- do.call(rbind, genebodyfpkms)
    row.names(genebodyfpkms) <- 1:nrow(genebodyfpkms)
    reslist[['GeneBodyFPKMstats']] <- genebodyfpkms
    
    reslist[['GeneBodyfwdFPMcum']] <- gbfwdcumfpm
    reslist[['GeneBodyrevFPMcum']] <- gbrevcumfpm
    
    reslist[['GeneBodysinglefwdFPM']] <- gbfwdindividual
    reslist[['GeneBodysinglerevFPM']] <- gbrevindividual
    
    
  }
  
  return(reslist)
  
}


combinesingle <- function(tssradius = tssradius, ttsradius = ttsradius,
                          singlegenename = 'B9d1',
                          TSSsingleFPMs = TSSsinglefwdFPMs, TTSsingleFPMs = TTSsinglefwdFPMs,
                          GeneBodysingleFPMs = GeneBodysinglefwdFPMs){
  
  if(max(tssradius) > 0){
    TSSsinglepart <- as.numeric(TSSsingleFPMs[[singlegenename]])[1:max(tssradius)]
  }else{
    TSSsinglepart <- NULL
  }
  
  if(max(ttsradius) > 0){
    TTSsinglepart <- as.numeric(TTSsingleFPMs[[singlegenename]])[(max(ttsradius) + 1 + 1):
                                                                   length(as.numeric(TTSsingleFPMs[[singlegenename]]))]
  }else{
    TTSsinglepart <- NULL
  }
  
  combinesingleFPMmeans <- c(TSSsinglepart, as.numeric(GeneBodysingleFPMs[[singlegenename]]),
                             TTSsinglepart)
  
  return(combinesingleFPMmeans)
  
  
}


plotends <- function(fwdreads = TSSfwdFPMmeans, revreads = TSSrevFPMmeans,
                     widthlimit = max(tssradius), halfwidth = min(tssradius), endtype = 'TSS',
                     label = label){
  
  droppart <- 0
  if(!is.null(halfwidth)){
    if(halfwidth < widthlimit){
      droppart <- widthlimit - halfwidth
    }
  }else{
    halfwidth <- widthlimit
    droppart <- 0
  }
  
  if(is.vector(fwdreads)){
    
    fwdreads <- data.frame(reads = fwdreads, stringsAsFactors = FALSE)
    
    
  }
  fwdparts <- fwdreads[(droppart + 1):(nrow(fwdreads) - droppart),,drop = FALSE]
  #drop For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension.
  row.names(fwdparts) <- 1:nrow(fwdparts)
  
  if(!is.null(revreads)){
    if(is.vector(revreads)){
      revreads <- data.frame(reads = revreads, stringsAsFactors = FALSE)
      
    }
    revparts <- revreads[(droppart + 1):(nrow(revreads) - droppart),,drop = FALSE]
    #drop For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension.
    row.names(revparts) <- 1:nrow(revparts)
  }
  
  
  midcoord <- halfwidth + 1
  endcoord <- 2*halfwidth + 1
  startcoord <- 1
  
  i <- 1
  for(i in 1:ncol(fwdparts)){
    fwdpart <- fwdparts[,i, drop = FALSE]
    fwdpart$condition <- names(fwdpart)[1]
    fwdpart$strand <- 'sense'
    fwdpart$xcoord <- as.numeric(row.names(fwdpart))
    
    if(!is.null(revparts)){
      revpart <- revparts[,unique(fwdpart$condition), drop = FALSE]
      revpart[,1] <- -revpart[,1]
      revpart$condition <- names(revpart)[1]
      revpart$strand <- 'antisense'
      revpart$xcoord <- as.numeric(row.names(revpart))
    }
    
    reads_total <- rbind(fwdpart, revpart)
    names(reads_total)[1] <- 'reads'
    
    if(i == 1){
      reads_totals <- reads_total
    }else{
      reads_totals <- rbind(reads_totals, reads_total)
    }
    
  }
  
  reads_totals$strand <- factor(reads_totals$strand, levels = c('sense', 'antisense'), ordered = TRUE)
  
  plottitle <- paste0('Range around ', toupper(endtype), ' ', halfwidth, 'bp')
  
  if(!is.null(label)){
    plottitle <- paste0('Range around', toupper(endtype), ' ', halfwidth, 'bp (', label, ')')
    
  }
  
  library(ggplot2)
  library(scales)
  library(grid)
  
  #grid.newpage()
  p <- ggplot(reads_totals, mapping = aes(x = xcoord, y = reads, color = condition))
  
  p <- p + geom_line(size = 1) +
    scale_x_continuous(breaks = c(startcoord, midcoord, endcoord),
                       labels = c(paste0('-', halfwidth), toupper(endtype),
                                  paste0('+', halfwidth))) +
    xlab('') + ylab('FPM') +
    geom_vline(xintercept = midcoord, linetype = 2, color = 'red', size = 1) +
    facet_grid(strand~., scales = 'free_y') +
    ggtitle(plottitle) +
    scale_color_discrete(name = 'condition') +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.x=element_text(hjust=1, angle = 90)) +
    theme(panel.spacing = unit(0, 'lines'))
  
  if(length(unique(reads_totals$condition)) == 1){
    p <- p + guides(color = FALSE)
    
  }
  
  print(p)
  
}

plotregions <- function(fwdreads = combinefwdFPMmeans, revreads = combinerevFPMmeans,
                        upwidthlimit = max(tssradius), upwidth = min(tssradius),
                        dnwidthlimit = max(ttsradius), dnwidth = min(ttsradius),
                        endtype = c('TSS', 'TTS'), genesym = NULL,
                        label = label, pausepoint = NULL){
  
  updroppart <- 0
  if(!is.null(upwidth)){
    if(upwidth < upwidthlimit){
      updroppart <- upwidthlimit - upwidth
    }
  }else{
    upwidth <- upwidthlimit
    updroppart <- 0
  }
  
  dndroppart <- 0
  if(!is.null(dnwidth)){
    if(dnwidth < dnwidthlimit){
      dndroppart <- dnwidthlimit - dnwidth
    }
  }else{
    dnwidth <- dnwidthlimit
    dndroppart <- 0
  }
  
  
  
  if(is.vector(fwdreads)){
    
    fwdreads <- data.frame(reads = fwdreads, stringsAsFactors = FALSE)
    
    
  }
  fwdparts <- fwdreads[(updroppart + 1):(nrow(fwdreads) - dndroppart),,drop = FALSE]
  #drop For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension.
  row.names(fwdparts) <- 1:nrow(fwdparts)
  
  if(!is.null(revreads)){
    if(is.vector(revreads)){
      revreads <- data.frame(reads = revreads, stringsAsFactors = FALSE)
      
    }
    revparts <- revreads[(updroppart + 1):(nrow(revreads) - dndroppart),,drop = FALSE]
    #drop For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension.
    row.names(revparts) <- 1:nrow(revparts)
  }
  
  uppoints <- c(1, upwidth + 1)
  dnpoints <- c((nrow(fwdparts) - dnwidth), nrow(fwdparts))
  
  if(uppoints[1] == uppoints[2]){
    uppointlabels <- rep(endtype[1], 2)
  }else{
    uppointlabels <- c(paste0('-', upwidth), endtype[1])
  }
  
  if(dnpoints[1] == dnpoints[2]){
    dnpointlabels <- rep(endtype[2], 2)
  }else{
    dnpointlabels <- c(endtype[2], paste0('+', dnwidth))
  }
  
  if(!is.null(pausepoint)){
    pausecoord <- uppoints[2] + pausepoint
    pausepointlabel <- paste0('+', pausepoint)
    
  }else{
    pausecoord <- NULL
    pausepointlabel <- NULL
  }
  
  points <- c(uppoints, pausecoord, dnpoints)
  pointlabels <- c(uppointlabels, pausepointlabel, dnpointlabels)
  
  
  i <- 1
  for(i in 1:ncol(fwdparts)){
    fwdpart <- fwdparts[,i, drop = FALSE]
    fwdpart$condition <- names(fwdpart)[1]
    fwdpart$strand <- 'sense'
    fwdpart$xcoord <- as.numeric(row.names(fwdpart))
    reads_total <- fwdpart
    
    if(!is.null(revreads)){
      revpart <- revparts[,unique(fwdpart$condition), drop = FALSE]
      revpart[,1] <- -revpart[,1]
      revpart$condition <- names(revpart)[1]
      revpart$strand <- 'antisense'
      revpart$xcoord <- as.numeric(row.names(revpart))
      reads_total <- rbind(fwdpart, revpart)
    }
    
    names(reads_total)[1] <- 'reads'
    
    if(i == 1){
      reads_totals <- reads_total
    }else{
      reads_totals <- rbind(reads_totals, reads_total)
    }
    
  }
  
  reads_totals$strand <- factor(reads_totals$strand, levels = c('sense', 'antisense'), ordered = TRUE)
  
  pointlabelslen <- length(pointlabels)
  
  plottitle <- paste0('Range from ', c(paste0(pointlabels[1], 'bp of '), '')
                      [c(points[1] != points[2], points[1] == points[2])],
                      pointlabels[2], ' to ',
                      c(paste0(pointlabels[pointlabelslen], 'bp of '), '')
                      [c(points[pointlabelslen - 1] != points[pointlabelslen],
                         points[pointlabelslen - 1] == points[pointlabelslen])],
                      pointlabels[pointlabelslen - 1])
  
  if(!is.null(genesym)){
    plottitle <- paste0(genesym, ' ', plottitle)
  }
  
  if(!is.null(label)){
    plottitle <- paste0(plottitle, ' (', label, ')')
  }
  
  
  library(ggplot2)
  library(scales)
  library(grid)
  
  #grid.newpage()
  p <- ggplot(reads_totals, mapping = aes(x = xcoord, y = reads, color = condition))
  
  p <- p + geom_line(size = 1) +
    scale_x_continuous(breaks = points,
                       labels = pointlabels) +
    xlab('') + ylab('FPM') +
    geom_vline(xintercept = points[2], linetype = 2, color = 'red', size = 1) +
    geom_vline(xintercept = points[3], linetype = 2, color = 'red', size = 1) +
    facet_grid(strand~., scales = 'free_y') +
    ggtitle(plottitle) +
    scale_color_discrete(name = 'condition') +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.x=element_text(hjust=1, angle = 90)) +
    theme(panel.spacing = unit(0, 'lines'))
  
  if(length(unique(reads_totals$condition)) == 1){
    p <- p + guides(color = FALSE)
    
  }
  
  if(!is.null(pausecoord)){
    p <- p +
      geom_vline(xintercept = pausecoord, linetype = 2, color = 'red', size = 1)
    
  }
  
  
  
  print(p)
  
  
}



metaplot <- function(metafile = metawt0file,
                     targetgenefile = metagenefile, genomename = NULL,
                     tssradius = c(2000, 1000, 500), ttsradius = c(2000, 1000, 500),
                     genebodylen = 4000,
                     label = NULL,
                     strandmethod = 1,
                     threads = 4,
                     savegenenames = c('Fyn', 'Lats1', 'Pfn1'), plotgenenames = TRUE, plotfig = TRUE,
                     genelencutoff = NULL,
                     fpkmcutoff = NULL){
  
  #Libraries required
  library(GenomicAlignments)
  
  #Prepare the genes to GRanges object#####
  if(!is.null(targetgenefile)){
    
    if(file.exists(targetgenefile)){
      
      genes <- read.table(targetgenefile, sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE, quote = '',
                          check.names = FALSE)
      
      genes$start <- genes$start + 1
      
      genes <- genes[(abs(genes$end - genes$start) + 1) > max(c(tssradius, ttsradius, genelencutoff)),]
      
      geneframe <- genes
      
      genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      genes <- GRanges(seqnames = genes$seqnames,
                       ranges = IRanges(start = genes$start, end = genes$end),
                       strand = genes$strand,
                       gene_id = genes$gene_id)
      genes <- sort(genes)
      
    }
    
    
  }else if(!is.null(genomename)){
    
    genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    genes <- genes[genes@ranges@width > max(c(tssradius, ttsradius, genelencutoff))]
    
    genes <- genes[genes$updis >= max(c(tssradius))]
    genes <- genes[genes$dndis >= max(c(ttsradius))]
    
    genes <- sort(genes)
    names(genes) <- 1:length(genes)
    genes <- keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
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
  
  #Read in the strand-specific paired-end data#####
  
  reads <- parsereadsfile(readsfile = metafile, strandmethod = strandmethod)
  readsdepth <- length(reads) #FRAGMENT counts
  
  #FPKM screen###
  if(!is.null(fpkmcutoff)){
    
    genes <- getfpkm(features = genes, reads = reads, readsdepth = readsdepth)
    genes <- genes[genes$fpkm > fpkmcutoff]
    genes <- sort(genes)
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
    
    starts <- subgeneranges@ranges@start - max(tssradius, ttsradius)
    starts[starts < 1] <- 1
    
    widthes <- subgeneranges@ranges@width + max(tssradius, ttsradius)
    ends <- starts + widthes - 1
    endlimits <- reads@seqinfo@seqlengths
    names(endlimits) <- reads@seqinfo@seqnames
    endlimits <- endlimits[as.character(subgeneranges@seqnames)]
    ends[ends > endlimits] <- endlimits[ends > endlimits]
    widthes <- ends - starts + 1
    
    subgeneranges <- GRanges(seqnames = subgeneranges@seqnames,
                             ranges = IRanges(start = starts, width = widthes), strand = '*')
    #Set strand as '*' here, so that both the reads on sense strand and antisense strand of a gene
    #can be selected next to plot metagene covering both sense and antisense strands
    
    subreads <- subsetByOverlaps(reads, subgeneranges)
    #Here, directly select reads from the paired fragment alignments, no need to separately from first reads
    #alignments and last reads alignments, and convert the strand symbols for one of them from artifical
    #symbols to true symbols, because if select from paired fragment alignments, the strand symbols of these
    #fragments are from the reads alignments file with the true strand symbols, and the one with artifical
    #symbols is disregarded automatically.
    
    subreads <- sort(subreads)
    
    subreadsranges <- range(subreads)
    subgenes <- subsetByOverlaps(subgenes, subreadsranges)
    subgenes <- sort(subgenes)
    
    if(length(subreads) == 0 | length(subgenes) == 0){
      next()
    }
    
    unitlist <- list(subgenes = subgenes, subreads = subreads)
    parallellist[[t]] <- unitlist
    
  }
  
  if(length(parallellist) == 0){
    return(NULL)
  }
  
  #Use parallel#####
  #date()
  #subreports <- lapply(X = parallellist, FUN = mapfunction,
  #                     depth = readsdepth,
  #                     tssradius = tssradius, ttsradius = ttsradius, genebodylen = genebodylen,
  #                     saveindividual = savegenenames)
  #date()
  
  
  
  library(doParallel)
  #threads <- 8
  
  cores <- detectCores()
  cl <- makeCluster(min(threads, cores))
  registerDoParallel(cl)
  
  date()
  subreports <- foreach(datalist = parallellist,
                        .export = ls(name = globalenv())) %dopar% mapfunction(datalist,
                                                                              depth = readsdepth,
                                                                              tssradius = tssradius, ttsradius = ttsradius,
                                                                              genebodylen = genebodylen,
                                                                              saveindividual = savegenenames)
  date()
  
  stopCluster(cl)
  
  #Integrate#####
  gene_ids <- data.frame(gene_id = genes$gene_id, stringsAsFactors = FALSE)
  TSSFPKMstatses <- gene_ids
  TTSFPKMstatses <- gene_ids
  GeneBodyFPKMstatses <- gene_ids
  
  TSSfwdFPMcums <- NULL
  TSSrevFPMcums <- NULL
  TTSfwdFPMcums <- NULL
  TTSrevFPMcums <- NULL
  GeneBodyfwdFPMcums <- NULL
  GeneBodyrevFPMcums <- NULL
  combinefwdFPMcums <- NULL
  combinerevFPMcums <- NULL
  
  TSSsinglefwdFPMs <- NULL
  TSSsinglerevFPMs <- NULL
  TTSsinglefwdFPMs <- NULL
  GeneBodysinglefwdFPMs <- NULL
  TTSsinglerevFPMs <- NULL
  GeneBodysinglerevFPMs <- NULL
  
  combinesinglefwdFPMs <- list()
  combinesinglerevFPMs <- list()
  
  i <- 1
  for(i in 1:length(subreports)){
    
    subreport <- subreports[[i]]
    
    if('TSSFPKMstats' %in% names(subreport)){
      
      TSSFPKMstats <- subreport$TSSFPKMstats
      TSSfwdFPMcum <- subreport$TSSfwdFPMcum
      TSSrevFPMcum <- subreport$TSSrevFPMcum
      
      if(i == 1){
        TSSFPKMstatses <- TSSFPKMstats
        TSSfwdFPMcums <- TSSfwdFPMcum
        TSSrevFPMcums <- TSSrevFPMcum
        
      }else{
        TSSFPKMstatses <- rbind(TSSFPKMstatses, TSSFPKMstats)
        TSSfwdFPMcums <- TSSfwdFPMcums + TSSfwdFPMcum
        TSSrevFPMcums <- TSSrevFPMcums + TSSrevFPMcum
      }
      
      if('TSSsinglefwdFPM' %in% names(subreport)){
        
        TSSsinglefwdFPM <- subreport$TSSsinglefwdFPM
        TSSsinglerevFPM <- subreport$TSSsinglerevFPM
        
        if(i == 1){
          TSSsinglefwdFPMs <- TSSsinglefwdFPM
          TSSsinglerevFPMs <- TSSsinglerevFPM
        }else{
          TSSsinglefwdFPMs <- c(TSSsinglefwdFPMs, TSSsinglefwdFPM)
          TSSsinglerevFPMs <- c(TSSsinglerevFPMs, TSSsinglerevFPM)
        }
      }
      
      
    }
    
    if('TTSFPKMstats' %in% names(subreport)){
      
      TTSFPKMstats <- subreport$TTSFPKMstats
      TTSfwdFPMcum <- subreport$TTSfwdFPMcum
      TTSrevFPMcum <- subreport$TTSrevFPMcum
      
      if(i == 1){
        TTSFPKMstatses <- TTSFPKMstats
        TTSfwdFPMcums <- TTSfwdFPMcum
        TTSrevFPMcums <- TTSrevFPMcum
      }else{
        TTSFPKMstatses <- rbind(TTSFPKMstatses, TTSFPKMstats)
        TTSfwdFPMcums <- TTSfwdFPMcums + TTSfwdFPMcum
        TTSrevFPMcums <- TTSrevFPMcums + TTSrevFPMcum
      }
      
      if('TTSsinglefwdFPM' %in% names(subreport)){
        
        TTSsinglefwdFPM <- subreport$TTSsinglefwdFPM
        TTSsinglerevFPM <- subreport$TTSsinglerevFPM
        if(i == 1){
          TTSsinglefwdFPMs <- TTSsinglefwdFPM
          TTSsinglerevFPMs <- TTSsinglerevFPM
        }else{
          TTSsinglefwdFPMs <- c(TTSsinglefwdFPMs, TTSsinglefwdFPM)
          TTSsinglerevFPMs <- c(TTSsinglerevFPMs, TTSsinglerevFPM)
        }
        
      }
      
    }
    
    if('GeneBodyFPKMstats' %in% names(subreport)){
      
      GeneBodyFPKMstats <- subreport$GeneBodyFPKMstats
      GeneBodyfwdFPMcum <- subreport$GeneBodyfwdFPMcum
      GeneBodyrevFPMcum <- subreport$GeneBodyrevFPMcum
      
      if(i == 1){
        GeneBodyFPKMstatses <- GeneBodyFPKMstats
        GeneBodyfwdFPMcums <- GeneBodyfwdFPMcum
        GeneBodyrevFPMcums <- GeneBodyrevFPMcum
      }else{
        GeneBodyFPKMstatses <- rbind(GeneBodyFPKMstatses, GeneBodyFPKMstats)
        GeneBodyfwdFPMcums <- GeneBodyfwdFPMcums + GeneBodyfwdFPMcum
        GeneBodyrevFPMcums <- GeneBodyrevFPMcums + GeneBodyrevFPMcum
        
      }
      
      
      if(max(tssradius) > 0){
        TSSfwdpart <- TSSfwdFPMcums[1:max(tssradius)]
        TSSrevpart <- TSSrevFPMcums[1:max(tssradius)]
      }else{
        TSSfwdpart <- NULL
        TSSrevpart <- NULL
      }
      
      if(max(ttsradius) > 0){
        TTSfwdpart <- TTSfwdFPMcums[(max(ttsradius) + 1 + 1):length(TTSfwdFPMcums)]
        TTSrevpart <- TTSrevFPMcums[(max(ttsradius) + 1 + 1):length(TTSfwdFPMcums)]
      }else{
        TTSfwdpart <- NULL
        TTSrevpart <- NULL
      }
      
      combinefwdFPMcums <- c(TSSfwdpart, GeneBodyfwdFPMcums, TTSfwdpart)
      combinerevFPMcums <- c(TSSrevpart, GeneBodyrevFPMcums, TTSrevpart)
      
      
      
      if('GeneBodysinglefwdFPM' %in% names(subreport)){
        
        GeneBodysinglefwdFPM <- subreport$GeneBodysinglefwdFPM
        GeneBodysinglerevFPM <- subreport$GeneBodysinglerevFPM
        
        if(i == 1){
          GeneBodysinglefwdFPMs <- GeneBodysinglefwdFPM
          GeneBodysinglerevFPMs <- GeneBodysinglerevFPM
        }else{
          GeneBodysinglefwdFPMs <- c(GeneBodysinglefwdFPMs, GeneBodysinglefwdFPM)
          GeneBodysinglerevFPMs <- c(GeneBodysinglerevFPMs, GeneBodysinglerevFPM)
        }
      }
      
    }
    
    
  }
  
  FPKMstatses <- cbind(gene_ids, TSSFPKMstatses[-1], TTSFPKMstatses[-1], GeneBodyFPKMstatses[-1])
  
  TSSfwdFPMmeans <- TSSfwdFPMcums/nrow(gene_ids)
  TSSrevFPMmeans <- TSSrevFPMcums/nrow(gene_ids)
  TTSfwdFPMmeans <- TTSfwdFPMcums/nrow(gene_ids)
  TTSrevFPMmeans <- TTSrevFPMcums/nrow(gene_ids)
  GeneBodyfwdFPMmeans <- GeneBodyfwdFPMcums/nrow(gene_ids)
  GeneBodyrevFPMmeans <- GeneBodyrevFPMcums/nrow(gene_ids)
  combinefwdFPMmeans <- combinefwdFPMcums/nrow(gene_ids)
  combinerevFPMmeans <- combinerevFPMcums/nrow(gene_ids)
  
  if(max(tssradius) > 0 & plotfig == TRUE){
    
    j <- 1
    for(j in 1:length(tssradius)){
      tssr <- tssradius[j]
      
      plotends(fwdreads = TSSfwdFPMmeans, revreads = TSSrevFPMmeans,
               widthlimit = max(tssradius), halfwidth = tssr, endtype = 'TSS', label = label)
      
      
    }
    
  }
  
  if(max(ttsradius) > 0 & plotfig == TRUE){
    
    for(j in 1:length(ttsradius)){
      ttsr <- ttsradius[j]
      plotends(fwdreads = TTSfwdFPMmeans, revreads = TTSrevFPMmeans,
               widthlimit = max(ttsradius), halfwidth = ttsr, endtype = 'TTS', label = label)
      
    }
    
  }
  
  if(genebodylen > 0 & plotfig == TRUE){
    
    plotregions(fwdreads = GeneBodyfwdFPMmeans, revreads = GeneBodyrevFPMmeans,
                upwidthlimit = 0, upwidth = 0, dnwidthlimit = 0, dnwidth = 0,
                endtype = c('TSS', 'TTS'), label = label)
    
    plotregions(fwdreads = combinefwdFPMmeans, revreads = combinerevFPMmeans,
                upwidthlimit = max(tssradius), upwidth = NULL,
                dnwidthlimit = max(ttsradius), dnwidth = NULL,
                endtype = c('TSS', 'TTS'), label = label)
    
  }
  
  if(!is.null(savegenenames)){
    
    j <- 1
    for(j in 1:length(savegenenames)){
      
      plotsinglegenename <- savegenenames[j]
      
      if(genebodylen > 0){
        singlegenefwdcombine <- combinesingle(tssradius = tssradius, ttsradius = ttsradius,
                                              singlegenename = plotsinglegenename,
                                              TSSsingleFPMs = TSSsinglefwdFPMs,
                                              TTSsingleFPMs = TTSsinglefwdFPMs,
                                              GeneBodysingleFPMs = GeneBodysinglefwdFPMs)
        
        singlegenerevcombine <- combinesingle(tssradius = tssradius, ttsradius = ttsradius,
                                              singlegenename = plotsinglegenename,
                                              TSSsingleFPMs = TSSsinglerevFPMs,
                                              TTSsingleFPMs = TTSsinglerevFPMs,
                                              GeneBodysingleFPMs = GeneBodysinglerevFPMs)
        
        combinesinglefwdFPMs[[plotsinglegenename]] <- Rle(singlegenefwdcombine)
        combinesinglerevFPMs[[plotsinglegenename]] <- Rle(singlegenerevcombine)
        
        
        if(plotgenenames == TRUE & plotfig == TRUE){
          
          print(
            
            plotregions(fwdreads = singlegenefwdcombine, revreads = singlegenerevcombine,
                        upwidthlimit = max(tssradius), upwidth = NULL,
                        dnwidthlimit = max(ttsradius), dnwidth = NULL,
                        endtype = c('TSS', 'TTS'), genesym = plotsinglegenename, label = label)
            
          )
          
        }
        
        
        
      }
      
      
    }
    
  }
  
  
  reslist <- list()
  reslist[['TSSFPKMstatses']] <- TSSFPKMstatses
  reslist[['TTSFPKMstatses']] <- TTSFPKMstatses
  reslist[['GeneBodyFPKMstatses']] <- GeneBodyFPKMstatses
  
  reslist[['TSSfwdFPMmeans']] <- TSSfwdFPMmeans
  reslist[['TSSrevFPMmeans']] <- TSSrevFPMmeans
  reslist[['TTSfwdFPMmeans']] <- TTSfwdFPMmeans
  reslist[['TTSrevFPMmeans']] <- TTSrevFPMmeans
  reslist[['GeneBodyfwdFPMmeans']] <- GeneBodyfwdFPMmeans
  reslist[['GeneBodyrevFPMmeans']] <- GeneBodyrevFPMmeans
  
  reslist[['combinefwdFPMmeans']] <- combinefwdFPMmeans
  reslist[['combinerevFPMmeans']] <- combinerevFPMmeans
  
  reslist[['TSSsinglefwdFPMs']] <- TSSsinglefwdFPMs
  reslist[['TSSsinglerevFPMs']] <- TSSsinglerevFPMs
  reslist[['TTSsinglefwdFPMs']] <- TTSsinglefwdFPMs
  reslist[['TTSsinglerevFPMs']] <- TTSsinglerevFPMs
  reslist[['GeneBodysinglefwdFPMs']] <- GeneBodysinglefwdFPMs
  reslist[['GeneBodysinglerevFPMs']] <- GeneBodysinglerevFPMs
  
  reslist[['combinesinglefwdFPMs']] <- combinesinglefwdFPMs
  reslist[['combinesinglerevFPMs']] <- combinesinglerevFPMs
  
  return(reslist)
  
}


mmetaplot <- function(metafiles = c(metawt0file, metako0file),
                      targetgenefile = metagenefile,
                      genomename = NULL,
                      tssradius = c(2000, 1000, 500), ttsradius = c(2000, 1000, 500),
                      genebodylen = 4000,
                      labels = c('WT', 'KO'),
                      strandmethod = 1,
                      threads = 4,
                      savegenenames = c('Fyn', 'Lats1', 'Pfn1'), plotgenenames = TRUE,
                      mergecases = FALSE,
                      genelencutoff = genelencutoff,
                      fpkmcutoff = NULL){
  
  
  metafiles <- getreadslist(timefiles = metafiles, mergefiles = mergecases, strandmode = strandmethod)
  metafilelength <- length(metafiles)
  
  reslist <- list()
  i <- 1
  for(i in 1:metafilelength){
    
    metafile <- metafiles[[i]]
    
    label <- labels[i]
    
    res <- metaplot(metafile = metafile,
                    targetgenefile = targetgenefile, genomename = genomename,
                    tssradius = tssradius, ttsradius = ttsradius,
                    genebodylen = genebodylen,
                    label = label,
                    strandmethod = strandmethod,
                    threads = threads,
                    savegenenames = savegenenames, plotgenenames = plotgenenames, plotfig = FALSE,
                    genelencutoff = genelencutoff, fpkmcutoff = fpkmcutoff)
    
    reslist[[label]] <- res
    
  }
  
  j <- 1
  for(j in 1:length(reslist)){
    res <- reslist[[j]]
    
    TSSFPKMstatses <- res$TSSFPKMstatses
    TTSFPKMstatses <- res$TTSFPKMstatses
    GeneBodyFPKMstatses <- res$GeneBodyFPKMstatses
    
    if(j == 1){
      mTSSfwdFPMmeans <- res$TSSfwdFPMmeans
      mTSSrevFPMmeans <- res$TSSrevFPMmeans
      mTTSfwdFPMmeans <- res$TTSfwdFPMmeans
      mTTSrevFPMmeans <- res$TTSrevFPMmeans
      mGeneBodyfwdFPMmeans <- res$GeneBodyfwdFPMmeans
      mGeneBodyrevFPMmeans <- res$GeneBodyrevFPMmeans
      mcombinefwdFPMmeans <- res$combinefwdFPMmeans
      mcombinerevFPMmeans <- res$combinerevFPMmeans
      
      mTSSFPKMstatses <- TSSFPKMstatses
      mTTSFPKMstatses <- TTSFPKMstatses
      mGeneBodyFPKMstatses <- GeneBodyFPKMstatses
      
    }else{
      mTSSfwdFPMmeans <- cbind(mTSSfwdFPMmeans, res$TSSfwdFPMmeans)
      mTSSrevFPMmeans <- cbind(mTSSrevFPMmeans, res$TSSrevFPMmeans)
      mTTSfwdFPMmeans <- cbind(mTTSfwdFPMmeans, res$TTSfwdFPMmeans)
      mTTSrevFPMmeans <- cbind(mTTSrevFPMmeans, res$TTSrevFPMmeans)
      mGeneBodyfwdFPMmeans <- cbind(mGeneBodyfwdFPMmeans, res$GeneBodyfwdFPMmeans)
      mGeneBodyrevFPMmeans <- cbind(mGeneBodyrevFPMmeans, res$GeneBodyrevFPMmeans)
      mcombinefwdFPMmeans <- cbind(mcombinefwdFPMmeans, res$combinefwdFPMmeans)
      mcombinerevFPMmeans <- cbind(mcombinerevFPMmeans, res$combinerevFPMmeans)
      
      mTSSFPKMstatses <- cbind(mTSSFPKMstatses, TSSFPKMstatses[-1])
      mTTSFPKMstatses <- cbind(mTTSFPKMstatses, TTSFPKMstatses[-1])
      mGeneBodyFPKMstatses <- cbind(mGeneBodyFPKMstatses, GeneBodyFPKMstatses[-1])
      
    }
    
    
  }
  
  
  if(max(tssradius) > 0){
    
    
    mTSSfwdFPMmeans <- as.data.frame(mTSSfwdFPMmeans, stringsAsFactors = FALSE)
    mTSSrevFPMmeans <- as.data.frame(mTSSrevFPMmeans, stringsAsFactors = FALSE)
    colnames(mTSSfwdFPMmeans) <- colnames(mTSSrevFPMmeans) <- names(reslist)
    rownames(mTSSfwdFPMmeans) <- rownames(mTSSrevFPMmeans) <- 1:nrow(mTSSfwdFPMmeans)
    
    j <- 1
    for(j in 1:length(tssradius)){
      tssr <- tssradius[j]
      
      plotends(fwdreads = mTSSfwdFPMmeans, revreads = mTSSrevFPMmeans,
               widthlimit = max(tssradius), halfwidth = tssr, endtype = 'TSS', label = NULL)
      
      
    }
    
  }
  
  if(max(ttsradius) > 0){
    
    
    mTTSfwdFPMmeans <- as.data.frame(mTTSfwdFPMmeans, stringsAsFactors = FALSE)
    mTTSrevFPMmeans <- as.data.frame(mTTSrevFPMmeans, stringsAsFactors = FALSE)
    colnames(mTTSfwdFPMmeans) <- colnames(mTTSrevFPMmeans) <- names(reslist)
    rownames(mTTSfwdFPMmeans) <- rownames(mTTSrevFPMmeans) <- 1:nrow(mTTSfwdFPMmeans)
    
    j <- 1
    for(j in 1:length(ttsradius)){
      ttsr <- ttsradius[j]
      
      plotends(fwdreads = mTTSfwdFPMmeans, revreads = mTTSrevFPMmeans,
               widthlimit = max(ttsradius), halfwidth = ttsr, endtype = 'TTS', label = NULL)
      
      
    }
    
  }
  
  if(genebodylen > 0){
    
    
    mGeneBodyfwdFPMmeans <- as.data.frame(mGeneBodyfwdFPMmeans, stringsAsFactors = FALSE)
    mGeneBodyrevFPMmeans <- as.data.frame(mGeneBodyrevFPMmeans, stringsAsFactors = FALSE)
    colnames(mGeneBodyfwdFPMmeans) <- colnames(mGeneBodyrevFPMmeans) <- names(reslist)
    rownames(mGeneBodyfwdFPMmeans) <- rownames(mGeneBodyrevFPMmeans) <- 1:nrow(mGeneBodyfwdFPMmeans)
    
    plotregions(fwdreads = mGeneBodyfwdFPMmeans, revreads = mGeneBodyrevFPMmeans,
                upwidthlimit = 0, upwidth = 0, dnwidthlimit = 0, dnwidth = 0,
                endtype = c('TSS', 'TTS'), label = label)
    
    
    
    mcombinefwdFPMmeans <- as.data.frame(mcombinefwdFPMmeans, stringsAsFactors = FALSE)
    mcombinerevFPMmeans <- as.data.frame(mcombinerevFPMmeans, stringsAsFactors = FALSE)
    colnames(mcombinefwdFPMmeans) <- colnames(mcombinerevFPMmeans) <- names(reslist)
    rownames(mcombinefwdFPMmeans) <- rownames(mcombinerevFPMmeans) <- 1:nrow(mcombinefwdFPMmeans)
    
    plotregions(fwdreads = mcombinefwdFPMmeans, revreads = mcombinerevFPMmeans,
                upwidthlimit = max(tssradius), upwidth = NULL,
                dnwidthlimit = max(ttsradius), dnwidth = NULL,
                endtype = c('TSS', 'TTS'), label = label)
    
  }
  
  
  
  if(!is.null(savegenenames) & plotgenenames == TRUE){
    
    k <- 1
    for(k in 1:length(savegenenames)){
      
      plotsinglegenename <- savegenenames[k]
      
      if(genebodylen > 0){
        
        l <- 1
        for(l in 1:length(reslist)){
          res <- reslist[[l]]
          
          if(l == 1){
            mcombinesinglefwdFPMs <- as.numeric(res$combinesinglefwdFPMs[[plotsinglegenename]])
            mcombinesinglerevFPMs <- as.numeric(res$combinesinglerevFPMs[[plotsinglegenename]])
          }else{
            mcombinesinglefwdFPMs <- cbind(mcombinesinglefwdFPMs,
                                           as.numeric(res$combinesinglefwdFPMs[[plotsinglegenename]]))
            mcombinesinglerevFPMs <- cbind(mcombinesinglerevFPMs,
                                           as.numeric(res$combinesinglerevFPMs[[plotsinglegenename]]))
          }
          
        }
        
        if(is.vector(mcombinesinglefwdFPMs)){
          mcombinesinglefwdFPMs <- matrix(mcombinesinglefwdFPMs, ncol = 1, byrow = FALSE)
          
        }
        
        if(is.vector(mcombinesinglerevFPMs)){
          mcombinesinglerevFPMs <- matrix(mcombinesinglerevFPMs, ncol = 1, byrow = FALSE)
        }
        
        
        colnames(mcombinesinglefwdFPMs) <- colnames(mcombinesinglerevFPMs) <- names(reslist)
        rownames(mcombinesinglefwdFPMs) <- rownames(mcombinesinglerevFPMs) <- 1:nrow(mcombinesinglefwdFPMs)
        mcombinesinglefwdFPMs <- as.data.frame(mcombinesinglefwdFPMs, stringsAsFactors = FALSE)
        mcombinesinglerevFPMs <- as.data.frame(mcombinesinglerevFPMs, stringsAsFactors = FALSE)
        
        if(ncol(mcombinesinglefwdFPMs) == 1){
          labelval <- colnames(mcombinesinglefwdFPMs)
        }else{
          labelval <- NULL
        }
        
        plotregions(fwdreads = mcombinesinglefwdFPMs, revreads = mcombinesinglerevFPMs,
                    upwidthlimit = max(tssradius), upwidth = NULL,
                    dnwidthlimit = max(ttsradius), dnwidth = NULL,
                    endtype = c('TSS', 'TTS'), genesym = plotsinglegenename, label = labelval)
        
        
      }
      
    }
    
    
  }
  
  if(!is.null(mTSSFPKMstatses)){
    
    names(mTSSFPKMstatses) <- c('gene_id', paste0(names(mTSSFPKMstatses)[-1],
                                                  '.',
                                                  rep(names(reslist),
                                                      each = ncol(reslist[[1]]$TSSFPKMstatses) - 1)))
    
  }
  
  if(!is.null(mTTSFPKMstatses)){
    
    names(mTTSFPKMstatses) <- c('gene_id', paste0(names(mTTSFPKMstatses)[-1],
                                                  '.',
                                                  rep(names(reslist),
                                                      each = ncol(reslist[[1]]$TTSFPKMstatses) - 1)))
    
  }
  
  if(!is.null(mGeneBodyFPKMstatses)){
    
    names(mGeneBodyFPKMstatses) <- c('gene_id', paste0(names(mGeneBodyFPKMstatses)[-1],
                                                       '.',
                                                       rep(names(reslist),
                                                           each = ncol(reslist[[1]]$GeneBodyFPKMstatses) - 1)))
  }
  
  reslist$TSSFPKMstatses <- mTSSFPKMstatses
  reslist$TTSFPKMstatses <- mTTSFPKMstatses
  reslist$GeneBodyFPKMstatses <- mGeneBodyFPKMstatses
  
  return(reslist)
  
}



#Pausing index#####

makepauseregion <- function(filename = pauseidxtargetfile, pauselen = 1000, genomename = 'mm10',
                            genelencutoff = NULL){
  
  #Prepare the genes to GRanges object#####
  if(!is.null(filename)){
    
    if(file.exists(filename)){
      
      genes <- read.table(filename, sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE, quote = '',
                          check.names = FALSE)
      
      genes$start <- genes$start + 1
      
      genes <- genes[(abs(genes$end - genes$start) + 1) > max(pauselen, genelencutoff),]
      
      geneframe <- genes
      
      genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      
      genes$tss <- NA
      genes$tts <- NA
      genes$pausestart <- NA
      genes$pauseend <- NA
      genes$elongstart <- NA
      genes$elongend <- NA
      for(i in 1:nrow(genes)){
        
        if(genes[i, 4] == '+'){
          genes[i, 6] <- genes[i, 2]
          genes[i, 7] <- genes[i, 3]
          genes[i, 8] <- genes[i, 6] - pauselen
          genes[i, 9] <- genes[i, 6] + pauselen
          genes[i, 10] <- genes[i, 6] + pauselen + 1
          genes[i, 11] <- genes[i, 7]
        }else{
          genes[i, 6] <- genes[i, 3]
          genes[i, 7] <- genes[i, 2]
          genes[i, 8] <- genes[i, 6] - pauselen
          genes[i, 9] <- genes[i, 6] + pauselen
          genes[i, 10] <- genes[i, 7]
          genes[i, 11] <- genes[i, 6] - pauselen - 1
        }
        
      }
      
      genes <- genes[order(genes$chr, genes$strand, genes$tss, genes$tts),]
      row.names(genes) <- 1:nrow(genes)
      
      pause <- genes[c('chr', 'pausestart', 'pauseend',
                       'strand', 'gene_id')]
      elong <- genes[c('chr', 'elongstart', 'elongend',
                       'strand', 'gene_id')]
      
      names(pause) <- names(elong) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      
      pause <- GRanges(seqnames = pause$seqnames,
                       ranges = IRanges(start = pause$start, end = pause$end),
                       strand = pause$strand,
                       gene_id = pause$gene_id)
      
      elong <- GRanges(seqnames = elong$seqnames,
                       ranges = IRanges(start = elong$start, end = elong$end),
                       strand = elong$strand,
                       gene_id = elong$gene_id)
      
      genebody <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genebody) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      genebody <- GRanges(seqnames = genebody$seqnames,
                          ranges = IRanges(start = genebody$start, end = genebody$end),
                          strand = genebody$strand,
                          gene_id = genebody$gene_id)
      
    }
    
    
  }else if(!is.null(genomename)){
    
    genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    genes <- genes[genes@ranges@width > max(pauselen, genelencutoff)]
    
    genes <- genes[genes$updis >= pauselen]
    
    genes <- sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)), pruning.mode = 'coarse')
    
    genes$updis <- NULL
    genes$dndis <- NULL
    
    geneframe <- as.data.frame(genes)
    
    geneframe <- geneframe[c('seqnames', 'start', 'end', 'strand', 'gene_id')]
    
    geneframe$tss <- NA
    geneframe$tts <- NA
    geneframe$pausestart <- NA
    geneframe$pauseend <- NA
    geneframe$elongstart <- NA
    geneframe$elongend <- NA
    
    for(i in 1:nrow(geneframe)){
      
      if(geneframe[i, 4] == '+'){
        geneframe[i, 6] <- geneframe[i, 2]
        geneframe[i, 7] <- geneframe[i, 3]
        geneframe[i, 8] <- geneframe[i, 6] - pauselen
        geneframe[i, 9] <- geneframe[i, 6] + pauselen
        geneframe[i, 10] <- geneframe[i, 6] + pauselen + 1
        geneframe[i, 11] <- geneframe[i, 7]
      }else{
        geneframe[i, 6] <- geneframe[i, 3]
        geneframe[i, 7] <- geneframe[i, 2]
        geneframe[i, 8] <- geneframe[i, 6] - pauselen
        geneframe[i, 9] <- geneframe[i, 6] + pauselen
        geneframe[i, 10] <- geneframe[i, 7]
        geneframe[i, 11] <- geneframe[i, 6] - pauselen - 1
      }
      
    }
    
    pause <- geneframe[c('seqnames', 'pausestart', 'pauseend',
                         'strand', 'gene_id')]
    elong <- geneframe[c('seqnames', 'elongstart', 'elongend',
                         'strand', 'gene_id')]
    
    names(pause) <- names(elong) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
    
    
    pause <- GRanges(seqnames = pause$seqnames,
                     ranges = IRanges(start = pause$start, end = pause$end),
                     strand = pause$strand,
                     gene_id = pause$gene_id)
    
    elong <- GRanges(seqnames = elong$seqnames,
                     ranges = IRanges(start = elong$start, end = elong$end),
                     strand = elong$strand,
                     gene_id = elong$gene_id)
    
    genebody <- geneframe[c('seqnames', 'start', 'end', 'strand', 'gene_id')]
    names(genebody) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
    genebody <- GRanges(seqnames = genebody$seqnames,
                        ranges = IRanges(start = genebody$start, end = genebody$end),
                        strand = genebody$strand,
                        gene_id = genebody$gene_id)
    
    
  }
  
  reslist <- list()
  reslist$pause <- pause
  reslist$elong <- elong
  reslist$genebody <- genebody
  
  return(reslist)
  
}

pausefunction <- function(datalist = parallellist[[3]], depth = readsdepth,
                          saveindividual = savegenenames){
  
  #Convert the bam file to strand separated wig files
  strandsign <- as.character(datalist$subgenebody@strand@values)
  
  if('+' %in% strandsign){
    
    Fp <- towig(reads = datalist$subreads, strand = '+')
    
  }
  
  if('-' %in% strandsign){
    
    Fm <- towig(reads = datalist$subreads, strand = '-')
    
  }
  
  if(!('+' %in% strandsign)){
    
    Fp <- NULL
    
  }
  
  if(!('-' %in% strandsign)){
    
    Fm <- NULL
    
  }
  
  
  syms <- c()
  dens1s <- c()
  dens2s <- c()
  idxes <- c()
  fishers <- c()
  
  pauselist <- list()
  elonglist <- list()
  combinelist <- list()
  
  pause <- datalist$subpasue
  elong <- datalist$subelong
  
  pause <- as.data.frame(pause, stringsAsFactors = FALSE)
  elong <- as.data.frame(elong, stringsAsFactors = FALSE)
  
  i <- 1
  for(i in 1:nrow(pause)){
    #print(i)
    
    sym <- pause[i, 6]
    start1 <- pause[i, 2]
    end1 <- pause[i, 3]
    start2 <- elong[i, 2]
    end2 <- elong[i, 3]
    len1 <- end1 - start1 + 1
    len2 <- end2 - start2 + 1
    
    if(pause[i, 5] == '+') {
      emis1  <-
        (as.numeric(Fp[[as.character(pause[i, 1])]]))[c(start1:end1)]
      emis2  <-
        (as.numeric(Fp[[as.character(elong[i, 1])]]))[c(start2:end2)]
    }else {
      emis1  <-
        rev((as.integer(Fm[[as.character(pause[i, 1])]]))
            [c(start1:end1)])
      #Note here the minus strand has been reversed!
      emis2  <-
        rev((as.integer(Fm[[as.character(elong[i, 1])]]))
            [c(start2:end2)])
      #Note here the minus strand has been reversed!
    }
    
    partdepth1 <- sum(emis1)
    partdepth2 <- sum(emis2)
    
    if(partdepth1 == 0 | partdepth2 == 0){
      next()
    }
    
    emis1_adj <- emis1/(depth/10^6)
    emis2_adj <- emis2/(depth/10^6)
    
    
    if(sym %in% saveindividual){
      pauselist[[sym]] <- Rle(emis1_adj)
      elonglist[[sym]] <- Rle(emis2_adj)
      combinelist[[sym]] <- Rle(c(emis1_adj, emis2_adj))
    }
    
    #Add pseudo count 1!!!
    emis1_adj <- emis1_adj + 1
    emis2_adj <- emis2_adj + 1
    
    
    depth1 <- sum(emis1_adj)
    depth2 <- sum(emis2_adj)
    
    dens1 <- depth1/len1
    dens2 <- depth2/len2
    pauseidx <- dens1/dens2
    
    expec1 <- round((depth1 + depth2)*len1/(len1 + len2))
    expec2 <- round((depth1 + depth2)*len2/(len1 + len2))
    fishermatrix <- matrix(c(round(depth1), round(depth2), expec1, expec2),
                           nrow = 2, byrow = FALSE)
    fisher <- fisher.test(fishermatrix)$p.value
    
    syms <- c(syms, sym)
    dens1s <- c(dens1s, dens1)
    dens2s <- c(dens2s, dens2)
    idxes <- c(idxes, pauseidx)
    fishers <- c(fishers, fisher)
    
  }
  
  report <- data.frame(gene_id = syms,
                       pausedens = dens1s,
                       elongdens = dens2s,
                       pauseidx = idxes,
                       pval = fishers,
                       stringsAsFactors = FALSE)
  
  reslist <- list()
  reslist$report <- report
  reslist$singlepause <- pauselist
  reslist$singleelong <- elonglist
  reslist$singlecombine <- combinelist
  
  return(reslist)
  
}

calpauseidx <- function(bamfile = time1chrsfile,
                        genefile = pauseidxtargetfile, genomename = NULL,
                        tssradius = 1000,
                        label = NULL, strandmethod = 1,
                        savegenenames = c('Nolc1', 'Lamp2', 'Pdk3'), plotgenenames = TRUE,
                        threads = 4,
                        genelencutoff = NULL,
                        fpkmcutoff = NULL){
  
  library(GenomicAlignments)
  
  #Prepare the genes to GRanges object
  regionlist <- makepauseregion(filename = genefile, pauselen = tssradius, genomename = genomename,
                                genelencutoff = genelencutoff)
  
  pause <- regionlist$pause
  elong <- regionlist$elong
  genebody <- regionlist$genebody
  
  #Read in the strand-specific paired-end data#####
  reads <- parsereadsfile(readsfile = bamfile, strandmethod = strandmethod)
  readsdepth <- length(reads) #FRAGMENT counts
  
  #FPKM screen###
  if(!is.null(fpkmcutoff)){
    
    genebody <- getfpkm(features = genebody, reads = reads, readsdepth = readsdepth)
    genebody <- genebody[genebody$fpkm > fpkmcutoff]
    genebody <- sort(genebody)
    names(genebody) <- 1:length(genebody)
    
    pause <- pause[match(genebody$gene_id, pause$gene_id)]
    names(pause) <- 1:length(pause)
    
    elong <- elong[match(genebody$gene_id, elong$gene_id)]
    names(elong) <- 1:length(elong)
    
  }
  
  
  #Prepare for parallel######
  if(threads > length(genebody)){
    threads <- length(genebody)
  }
  
  subsize <- ceiling(length(genebody)/threads)
  
  parallellist <- list()
  t <- 4
  
  for(t in 1:threads){
    subpause <- pause[((t - 1)*subsize + 1):min(t*subsize, length(pause)),]
    subelong <- elong[((t - 1)*subsize + 1):min(t*subsize, length(elong)),]
    subgenebody <- genebody[((t - 1)*subsize + 1):min(t*subsize, length(genebody)),]
    
    subpauseranges <- range(subpause)
    subgenebodyranges <- range(subgenebody)
    
    
    
    subpausestarts <- subpauseranges@ranges@start
    substarts <- subgenebodyranges@ranges@start
    
    subpausewidthes <- subpauseranges@ranges@width
    subwidths <- subgenebodyranges@ranges@width
    
    subpauseends <- subpausestarts + subpausewidthes - 1
    subends <- substarts + subwidths - 1
    
    starts <- pmin(subpausestarts, substarts)
    ends <- pmax(subpauseends, subends)
    
    widthes <- ends - starts + 1
    
    subranges <- GRanges(seqnames = subgenebodyranges@seqnames,
                         ranges = IRanges(start = starts, width = widthes),
                         strand = subgenebodyranges@strand)
    
    subreads <- subsetByOverlaps(reads, subranges)
    
    subreads <- sort(subreads)
    
    subreadsranges <- range(subreads)
    
    
    subpause <- subsetByOverlaps(subpause, subreadsranges)
    subelong <- subsetByOverlaps(subelong, subreadsranges)
    subgenebody <- subsetByOverlaps(subgenebody, subreadsranges)
    
    if(length(subreads) == 0 | length(subpause) == 0){
      next()
    }
    
    unitlist <- list(subpasue = subpause, subelong = subelong, subgenebody = subgenebody,
                     subreads = subreads)
    parallellist[[t]] <- unitlist
    
  }
  
  if(length(parallellist) == 0){
    return(NULL)
  }
  
  #Use parallel#####
  #date()
  #subreports <- lapply(X = parallellist, FUN = pausefunction,
  #                     depth = readsdepth,
  #                     saveindividual = savegenenames)
  #date()
  
  library(doParallel)
  #threads <- 1
  
  cores <- detectCores()
  cl <- makeCluster(min(threads, cores))
  registerDoParallel(cl)
  
  date()
  subreports <- foreach(datalist = parallellist,
                        .export = ls(name = globalenv())) %dopar% pausefunction(datalist,
                                                                                depth = readsdepth,
                                                                                saveindividual = savegenenames)
  date()
  
  stopCluster(cl)
  
  
  
  
  
  for(i in 1:length(subreports)){
    if(i == 1){
      report <- subreports[[i]]$report
      singlepause <- subreports[[i]]$singlepause
      singleelong <- subreports[[i]]$singleelong
      singlecombine <- subreports[[i]]$singlecombine
      
    }else{
      report <- rbind(report, subreports[[i]]$report)
      singlepause <- c(singlepause, subreports[[i]]$singlepause)
      singleelong <- c(singleelong, subreports[[i]]$singleelong)
      singlecombine <- c(singlecombine, subreports[[i]]$singlecombine)
    }
  }
  
  report$padj <- p.adjust(report$pval, method = 'BH')
  
  
  if(!is.null(savegenenames) & plotgenenames == TRUE){
    
    j <- 1
    for(j in 1:length(savegenenames)){
      
      plotsinglegenename <- savegenenames[j]
      
      
      print(
        
        plotregions(fwdreads = as.numeric(singlecombine[[plotsinglegenename]]),
                    revreads = NULL,
                    upwidthlimit = tssradius, upwidth = NULL,
                    dnwidthlimit = 0, dnwidth = NULL,
                    endtype = c('TSS', 'TTS'), genesym = plotsinglegenename, label = label,
                    pausepoint = tssradius)
        
      )
      
    }
    
    
  }
  
  res <- list()
  res$report <- report
  res$singlepause <- singlepause
  res$singleelong <- singleelong
  res$singlecombine <- singlecombine
  
  return(res)
  
}



mcalpauseidx <- function(bamfiles = c(time1chrsfile, time2chrsfile),
                         genefile = pauseidxtargetfile, genomename = NULL,
                         tssradius = 1000,
                         labels = c('0min', '15min'), strandmethod = 1,
                         savegenenames = c('Nolc1', 'Lamp2', 'Pdk3'), plotgenenames = TRUE,
                         mergecases = FALSE,
                         threads = 4,
                         genelencutoff = NULL,
                         fpkmcutoff = NULL){
  
  
  bamfiles <- getreadslist(timefiles = bamfiles, mergefiles = mergecases, strandmode = strandmethod)
  bamfilelength <- length(bamfiles)
  
  reslist <- list()
  i <- 1
  for(i in 1:bamfilelength){
    
    bamfile <- bamfiles[[i]]
    
    label <- labels[i]
    
    res <- calpauseidx(bamfile = bamfile,
                       genefile = genefile, genomename = genomename,
                       tssradius = tssradius,
                       label = label, strandmethod = strandmethod,
                       savegenenames = savegenenames, plotgenenames = plotgenenames,
                       threads = threads,
                       genelencutoff = genelencutoff,
                       fpkmcutoff = fpkmcutoff)
    
    reslist[[label]] <- res
    
  }
  
  
  if(!is.null(savegenenames) & plotgenenames == TRUE){
    
    k <- 1
    for(k in 1:length(savegenenames)){
      
      plotsinglegenename <- savegenenames[k]
      
      l <- 1
      for(l in 1:length(reslist)){
        res <- reslist[[l]]
        
        if(l == 1){
          msinglecombines <- as.numeric(res$singlecombine[[plotsinglegenename]])
        }else{
          msinglecombines <- cbind(msinglecombines,
                                   as.numeric(res$singlecombine[[plotsinglegenename]]))
        }
        
      }
      
      if(is.vector(msinglecombines)){
        msinglecombines <- matrix(msinglecombines, ncol = 1, byrow = FALSE)
        
      }
      
      
      colnames(msinglecombines) <- names(reslist)
      rownames(msinglecombines) <- 1:nrow(msinglecombines)
      msinglecombines <- as.data.frame(msinglecombines, stringsAsFactors = FALSE)
      
      if(ncol(msinglecombines) == 1){
        labelval <- colnames(msinglecombines)
      }else{
        labelval <- NULL
      }
      
      plotregions(fwdreads = msinglecombines, revreads = NULL,
                  upwidthlimit = tssradius, upwidth = NULL,
                  dnwidthlimit = 0, dnwidth = NULL,
                  endtype = c('TSS', 'TTS'), genesym = plotsinglegenename, label = labelval,
                  pausepoint = tssradius)
      
      
      
      
    }
    
    
  }
  
  return(reslist)
  
}



























#Gene structure#####

getgc <- function(targetfile, genomename = 'mm10', genenames = NULL){
  
  library(GenomicAlignments)
  
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
      
      genes <- GRanges(seqnames = genes$seqnames,
                       ranges = IRanges(start = genes$start, end = genes$end),
                       strand = genes$strand,
                       gene_id = genes$gene_id)
      
      genes <- genegccont(origenes = genes, genomename = genomename)
      
      genes <- sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    if(!is.null(genenames)){
      genes <- subset(genes, gene_id %in% genenames)
    }
    
    genes <- sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
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
  
  featureranges <- IRanges(start = featurestarts, end = featureends)
  featuregrs <- GRanges(seqnames = genes@seqnames, ranges = featureranges,
                        strand = strands, gene_id = genes$gene_id)
  
  return(featuregrs)
  
}

genekmer <- function(genes = genes, genomename = 'mm10', k = 6){
  
  library(GenomicRanges)
  library(Biostrings)
  
  geneends <- genes@ranges@start + genes@ranges@width - 1
  
  if(genomename == 'mm10'){
    library("BSgenome.Mmusculus.UCSC.mm10")
    chrlimits <- Mmusculus@seqinfo@seqlengths
    names(chrlimits) <- Mmusculus@seqinfo@seqnames
    
    genelimits <- as.vector(chrlimits[as.character(genes@seqnames)])
    geneends[geneends > genelimits] <- genelimits[geneends > genelimits]
    genes <- GRanges(seqnames = genes@seqnames,
                     ranges = IRanges(start = genes@ranges@start,
                                      end = geneends),
                     strand = genes@strand,
                     gene_id = genes$gene_id)
    
    seqs <- getSeq(Mmusculus, genes)
  }else if(genomename == 'hg38'){
    library("BSgenome.Hsapiens.UCSC.hg38")
    chrlimits <- Hsapiens@seqinfo@seqlengths
    names(chrlimits) <- Hsapiens@seqinfo@seqnames
    
    genelimits <- as.vector(chrlimits[as.character(genes@seqnames)])
    geneends[geneends > genelimits] <- genelimits[geneends > genelimits]
    genes <- GRanges(seqnames = genes@seqnames,
                     ranges = IRanges(start = genes@ranges@start,
                                      end = geneends),
                     strand = genes@strand,
                     gene_id = genes$gene_id)
    
    seqs <- getSeq(Hsapiens, genes)
  }
  
  counts_genes <- oligonucleotideFrequency(seqs, width = k, step = 1)
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

getkmer <- function(targetfile, genomename = 'mm10', k = 6,
                    genes1 = subset(intracompareres$absolute_rate, absq == 'quantile4')$gene_id,
                    genes2 = subset(intracompareres$absolute_rate, absq == 'quantile1')$gene_id,
                    feature = c('promoter', 'end', 'genebody'),
                    radius = 1000){
  
  library(GenomicAlignments)
  
  if(!is.null(targetfile)){
    
    if(file.exists(targetfile)){
      
      genes <- read.table(targetfile, sep = '\t', header = TRUE,
                          stringsAsFactors = FALSE, quote = '',
                          check.names = FALSE)
      
      genes$start <- genes$start + 1
      
      geneframe <- genes
      
      genes <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genes) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      
      genes <- GRanges(seqnames = genes$seqnames,
                       ranges = IRanges(start = genes$start, end = genes$end),
                       strand = genes$strand,
                       gene_id = genes$gene_id)
      
      genes <- sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    genes <- sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
                           pruning.mode = 'coarse')
    
    genes$updis <- NULL
    genes$dndis <- NULL
    genes$GC <- NULL
    genes$exon <- NULL
    
  }
  
  genes1 <- subset(genes, gene_id %in% genes1)
  genes2 <- subset(genes, gene_id %in% genes2)
  
  genes1 <- sort(genes1)
  genes2 <- sort(genes2)
  names(genes1) <- 1:length(genes1)
  names(genes2) <- 1:length(genes2)
  
  if(feature != 'genebody'){
    genes1 <- extractfeature(genes = genes1, feature = feature, radius = radius)
    genes2 <- extractfeature(genes = genes2, feature = feature, radius = radius)
  }else{
    genes1 <- genes1
    genes2 <- genes2
  }
  
  
  kmerres1 <- genekmer(genes = genes1, genomename = genomename, k = k)
  kmerres2 <- genekmer(genes = genes2, genomename = genomename, k = k)
  
  kmerres <- selectfisher(freqvec1 = kmerres2, freqvec2 = kmerres1)
  
  return(kmerres)
  
}



getexon <- function(targetfile, genomename = 'mm10', genenames = NULL){
  
  library(GenomicAlignments)
  
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
      
      genes <- GRanges(seqnames = genes$seqnames,
                       ranges = IRanges(start = genes$start, end = genes$end),
                       strand = genes$strand,
                       gene_id = genes$gene_id)
      
      genes <- geneexonfreq(origenes = genes, genomename = genomename)
      
      genes <- sort(genes)
      
    }
    
  }else if(!is.null(genomename)){
    
    genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    if(!is.null(genenames)){
      genes <- subset(genes, gene_id %in% genenames)
    }
    
    genes <- sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)),
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

