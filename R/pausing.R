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
      
      
      pause <- GenomicRanges::GRanges(seqnames = pause$seqnames,
                                      ranges = IRanges::IRanges(start = pause$start, end = pause$end),
                                      strand = pause$strand,
                                      gene_id = pause$gene_id)
      
      elong <- GenomicRanges::GRanges(seqnames = elong$seqnames,
                                      ranges = IRanges::IRanges(start = elong$start, end = elong$end),
                                      strand = elong$strand,
                                      gene_id = elong$gene_id)
      
      genebody <- genes[c('chr', 'start', 'end', 'strand', 'gene_id')]
      names(genebody) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
      genebody <- GenomicRanges::GRanges(seqnames = genebody$seqnames,
                                         ranges = IRanges::IRanges(start = genebody$start, end = genebody$end),
                                         strand = genebody$strand,
                                         gene_id = genebody$gene_id)
      
    }
    
    
  }else if(!is.null(genomename)){
    
    genes <- get(paste0(genomename, '.packagegenes'))
    #genes <- readRDS(paste0(genomename, '.packagegenes.rds'))
    
    genes <- GenomicRanges::GRanges(genes)
    
    genes <- genes[genes@ranges@width > max(pauselen, genelencutoff)]
    
    genes <- genes[genes$updis >= pauselen]
    
    genes <- GenomicAlignments::sort(genes)
    names(genes) <- 1:length(genes)
    
    genes <- GenomeInfoDb::keepSeqlevels(x = genes, value = unique(as.character(genes@seqnames)), pruning.mode = 'coarse')
    
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
    
    
    pause <- GenomicRanges::GRanges(seqnames = pause$seqnames,
                                    ranges = IRanges::IRanges(start = pause$start, end = pause$end),
                                    strand = pause$strand,
                                    gene_id = pause$gene_id)
    
    elong <- GenomicRanges::GRanges(seqnames = elong$seqnames,
                                    ranges = IRanges::IRanges(start = elong$start, end = elong$end),
                                    strand = elong$strand,
                                    gene_id = elong$gene_id)
    
    genebody <- geneframe[c('seqnames', 'start', 'end', 'strand', 'gene_id')]
    names(genebody) <- c('seqnames', 'start', 'end', 'strand', 'gene_id')
    genebody <- GenomicRanges::GRanges(seqnames = genebody$seqnames,
                                       ranges = IRanges::IRanges(start = genebody$start, end = genebody$end),
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
      pauselist[[sym]] <- S4Vectors::Rle(emis1_adj)
      elonglist[[sym]] <- S4Vectors::Rle(emis2_adj)
      combinelist[[sym]] <- S4Vectors::Rle(c(emis1_adj, emis2_adj))
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



#'Calculate Pause indeces for genes from a bam file
#'
#'Calculate Pause indeces for genes from a bam file (Pause index is the ratio 
#'of transcription Polymerase II signal density near a gene promoter to signal 
#'density within the gene body, such that higher Pause indices reflect a 
#'greater enrichment of promoter-paused Polymerase II), and also plot the FPM 
#'value of each bp for specific genes. 
#'
#'@param bamfile The bam file need to calculate the Pause indeces for genes. 
#'  Can be a string indicating the directory of the file, or a GAlignmentPairs 
#'  object, or a GRanges object from the original bam file. 
#'@param genefile The genes whose Pause indeces need to be calculated. If it 
#'  is NULL, all the genes in the genome specified by the parameter 
#'  \code{genomename} will be analyzed. If provied by the user, columns named 
#'  as chr, start, end, strand, and gene_id are required. All genes should be 
#'  longer than the maximum value of the parameters \code{tssradius} and 
#'  \code{genelencutoff}.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{genefile} is NULL.
#'@param tssradius A numeric value indicating the radius of the promoter 
#'  region centering around the TSS 
#'@param label A string that will be included in the titles of the specific 
#'  gene plots to indicate the experimental condition. 
#'@param strandmethod The strand specific method used when preparing the 
#'  sequencing library, can be 1 for directional ligation method and 2 for 
#'  dUTP method. If the sample is sequenced using a single strand method, set 
#'  it as 0.
#'@param savegenenames For which genes their concrete FPM value for each bp 
#'  position need to be saved, or plotted.
#'@param plotgenenames Whether to plot the FPM value for each bp position for 
#'  the genes provided by the parameter \code{savegenenames}.
#'@param threads Number of threads to do the parallelization. Default is 1.
#'@param genelencutoff The cutoff on gene length (bp). The default value is 
#'  NULL, but if it is set, only genes longer than this cutoff will be 
#'  considerred.
#'@param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM 
#'  value greater than the cutoff will be considerred. Default is 1. 
#'@return Will generate a list containing a data.frame recording the Pause 
#'  indeces of the genes, as well as plot the FPM value of each bp for the 
#'  specific genes indicated by the parameter \code{savegenenames}.
#'@export
calpauseidx <- function(bamfile,
                        genefile, genomename = NULL,
                        tssradius = 1000,
                        label = NULL, strandmethod = 1,
                        savegenenames = NULL, plotgenenames = TRUE,
                        threads = 1,
                        genelencutoff = NULL,
                        fpkmcutoff = 1){
  
  #library(GenomicAlignments)
  
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
    if(length(genebody) == 0){
      return(NULL)
    }
    genebody <- GenomicAlignments::sort(genebody)
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
    
    
    
    if((t - 1)*subsize + 1 > length(pause)){
      
      next()
      
    }
    
    
    
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
    
    subranges <- GenomicRanges::GRanges(seqnames = subgenebodyranges@seqnames,
                                        ranges = IRanges::IRanges(start = starts, width = widthes),
                                        strand = subgenebodyranges@strand)
    
    subreads <- IRanges::subsetByOverlaps(reads, subranges)
    
    subreads <- GenomicAlignments::sort(subreads)
    
    subreadsranges <- range(subreads)
    
    
    subpause <- IRanges::subsetByOverlaps(subpause, subreadsranges)
    subelong <- IRanges::subsetByOverlaps(subelong, subreadsranges)
    subgenebody <- IRanges::subsetByOverlaps(subgenebody, subreadsranges)
    
    if(length(subreads) == 0 | length(subpause) == 0){
      next()
    }
    
    
    
    #sharedgenes <- intersect(intersect(subpause$gene_id, subelong$gene_id), 
    #                         subgenebody$gene_id)
    
    #subpause <- subpause[subpause$gene_id %in% sharedgenes]
    #subelong <- subelong[subelong$gene_id %in% sharedgenes]
    #subgenebody <- subgenebody[subgenebody$gene_id %in% sharedgenes]
    
    
    
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
  
  if(threads == 1){
    
    subreports <- list()
    
    subreports[[1]] <- pausefunction(datalist = parallellist[[1]], 
                                     depth = readsdepth, 
                                     saveindividual = savegenenames)
  }else{
    
    #library(doParallel)
    #threads <- 1
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    doParallel::registerDoParallel(cl)
    
    date()
    `%dopar%` <- foreach::`%dopar%`
    subreports <- foreach::foreach(datalist = parallellist,
                                   .export = ls(name = globalenv())) %dopar% pausefunction(datalist,
                                                                                           depth = readsdepth,
                                                                                           saveindividual = savegenenames)
    date()
    
    parallel::stopCluster(cl)
    
    unregister_dopar()
    
  }
  
  
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



#'Calculate Pause indeces for genes for multiple bam files
#'
#'Calculate Pause indeces for genes for mulitple bam files (Pause index is the 
#'ratio of transcription Polymerase II signal density near a gene promoter to 
#'signal density within the gene body, such that higher Pause indices reflect 
#'a greater enrichment of promoter-paused Polymerase II), and also plot the 
#'FPM value of each bp for specific genes. 
#'
#'@param bamfiles The bam files need to calculate the Pause indeces for genes. 
#'  Should be a vector with elements as strings indicating the directores of 
#'  the files.
#'@param genefile The genes whose Pause indeces need to be calculated. If it 
#'  is NULL, all the genes in the genome specified by the parameter 
#'  \code{genomename} will be analyzed. If provied by the user, columns named 
#'  as chr, start, end, strand, and gene_id are required. All genes should be 
#'  longer than the maximum value of the parameters \code{tssradius} and 
#'  \code{genelencutoff}.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{genefile} is NULL.
#'@param tssradius A numeric value indicating the radius of the promoter 
#'  region centering around the TSS. 
#'@param labels A vector with elements as strings to be included in the titles 
#'  of the specific gene plots to indicate the experimental conditions of the 
#'  bam files. There should be no replicated elements in this vector.
#'@param strandmethod The strand specific method used when preparing the 
#'  sequencing library, can be 1 for directional ligation method and 2 for 
#'  dUTP method. If the sample is sequenced using a single strand method, set 
#'  it as 0.
#'@param savegenenames For which genes their concrete FPM value for each bp 
#'  position need to be saved, or plotted.
#'@param plotgenenames Whether to plot the FPM value for each bp position for 
#'  the genes provided by the parameter \code{savegenenames}.
#'@param mergecases Whether merge the data of all bam files together to one, 
#'  and then use it to calculate gene Pause indeces and plot genes. Default is 
#'  FALSE.
#'@param threads Number of threads to do the parallelization. Default is 1.
#'@param genelencutoff The cutoff on gene length (bp). The default value is 
#'  NULL, but if it is set, only genes longer than this cutoff will be 
#'  considerred.
#'@param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM 
#'  value greater than the cutoff will be considerred. Default is 1. 
#'@return Will generate a list with several sub-lists recording the Pause 
#'  indeces of the genes, as well as plot the FPM value of each bp for the 
#'  specific genes indicated by the parameter \code{savegenenames}.
#'@export
mcalpauseidx <- function(bamfiles,
                         genefile, genomename = NULL,
                         tssradius = 1000,
                         labels, strandmethod = 1,
                         savegenenames = NULL, plotgenenames = TRUE,
                         mergecases = FALSE,
                         threads = 1,
                         genelencutoff = NULL,
                         fpkmcutoff = 1){
  
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
      if(length(reslist) == 0){
        return(NULL)
      }
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


