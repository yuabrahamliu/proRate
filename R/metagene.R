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
  
  
  BiocGenerics::strand(reads1) <- firststrands
  BiocGenerics::strand(reads2) <- secondstrands
  #Use strand(reads1) to change the strand symbols, don't use reads1@strand, which will not work,
  #or rewrite reads1 <- GRanges(seqnames = reads1@seqnames, ranges = reads1@ranges,
  #strand = firststrands), which will lost other important information in the original reads1, such as
  #reads1@seqinfo, which includes chromosome length information
  
  reads <- c(reads1, reads2)
  reads <- GenomicAlignments::sort(reads)
  
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
        
        tssfwdindividual[[genecoord$gene_id]] <- S4Vectors::Rle(tssfwdfpm)
        tssrevindividual[[genecoord$gene_id]] <- S4Vectors::Rle(tssrevfpm)
        
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
        
        ttsfwdindividual[[genecoord$gene_id]] <- S4Vectors::Rle(ttsfwdfpm)
        ttsrevindividual[[genecoord$gene_id]] <- S4Vectors::Rle(ttsrevfpm)
        
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
        
        gbfwdindividual[[genecoord$gene_id]] <- S4Vectors::Rle(gbfwdfpm)
        gbrevindividual[[genecoord$gene_id]] <- S4Vectors::Rle(gbrevfpm)
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
  
  #library(ggplot2)
  #library(scales)
  #library(grid)
  
  #grid.newpage()
  p <- ggplot2::ggplot(reads_totals, mapping = ggplot2::aes(x = xcoord, y = reads, color = condition))
  
  p <- p + ggplot2::geom_Line(size = 1) +
    ggplot2::scale_x_continuous(breaks = c(startcoord, midcoord, endcoord),
                                labels = c(paste0('-', halfwidth), toupper(endtype),
                                           paste0('+', halfwidth))) +
    ggplot2::xlab('') + ggplot2::ylab('FPM') +
    ggplot2::geom_vline(xintercept = midcoord, linetype = 2, color = 'red', size = 1) +
    ggplot2::facet_grid(strand~., scales = 'free_y') +
    ggplot2::ggtitle(plottitle) +
    ggplot2::scale_color_discrete(name = 'condition') +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust=1, angle = 90)) +
    ggplot2::theme(panel.spacing = grid::unit(0, 'lines'))
  
  if(length(unique(reads_totals$condition)) == 1){
    p <- p + ggplot2::guides(color = FALSE)
    
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
  
  
  #library(ggplot2)
  #library(scales)
  #library(grid)
  
  #grid.newpage()
  p <- ggplot2::ggplot(reads_totals, mapping = ggplot2::aes(x = xcoord, y = reads, color = condition))
  
  p <- p + ggplot2::geom_Line(size = 1) +
    ggplot2::scale_x_continuous(breaks = points,
                                labels = pointlabels) +
    ggplot2::xlab('') + ggplot2::ylab('FPM') +
    ggplot2::geom_vline(xintercept = points[2], linetype = 2, color = 'red', size = 1) +
    ggplot2::geom_vline(xintercept = points[3], linetype = 2, color = 'red', size = 1) +
    ggplot2::facet_grid(strand~., scales = 'free_y') +
    ggplot2::ggtitle(plottitle) +
    ggplot2::scale_color_discrete(name = 'condition') +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(hjust=1, angle = 90)) +
    ggplot2::theme(panel.spacing = grid::unit(0, 'lines'))
  
  if(length(unique(reads_totals$condition)) == 1){
    p <- p + ggplot2::guides(color = FALSE)
    
  }
  
  if(!is.null(pausecoord)){
    p <- p +
      ggplot2::geom_vline(xintercept = pausecoord, linetype = 2, color = 'red', size = 1)
    
  }
  
  
  
  print(p)
  
  
}



#'Generate metagene plots from a bam file
#'
#'Generate metagene plots from a bam file, as well as other information 
#'including the FPM value of each bp position in the metagene plot, FPKM value 
#'in specific TSS, TTS, or gene body regions of each gene, etc. 
#'
#'@param metafile The bam file need to generate metagene plots. Can be a 
#'  string indicating the directory of the file, or a GAlignmentPairs object, 
#'  or a GRanges object from the original bam file. 
#'@param targetgenefile The genes whose FPM values need to be merged together 
#'  to generate a metagene plot. If it is NULL, all the genes in the genome 
#'  specified by the parameter \code{genomename} will be analyzed. If provied 
#'  by the user, columns named as chr, start, end, strand, and gene_id are 
#'  required. All genes should be longer than the maximum value of the 
#'  parameters \code{tssradius}, \code{ttsradius}, and \code{genelencutoff}.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetgenefile} is NULL.
#'@param tssradius If want to plot the metagene in TSS region, set this 
#'  parameter as a vector with each element corresponding to a radius of the 
#'  TSS region centering around the TSS site. Then, for each radius value, a 
#'  metagene plot will be generated.
#'@param ttsradius If want to plot the metagene in TTS region, set this 
#'  parameter as a vector with each element corresponding to a radius of the 
#'  TTS region centering around the TTS site. Then, for each radius value, a 
#'  metagene plot will be generated. 
#'@param genebodylen If want to plot the metagene in gene body region, set 
#'  this parameter as a specific numeric value (bp). Because different genes 
#'  have different gene body length, before merge their FPM values to a 
#'  metagene, their gene body length will be scaled to a uniform length, which 
#'  is set by this parameter \code{genebodylen}. Genes with a longer length 
#'  will be compressed to this gene length, while the ones with a shorter 
#'  length will be exteneded. 
#'@param label A string that will be included in the title of the metagene 
#'  plot to indicate the experimental condition. 
#'@param strandmethod The strand specific method used when preparing the 
#'  sequencing library, can be 1 for directional ligation method and 2 for 
#'  dUTP method. If the sample is sequenced using a single strand method, set 
#'  it as 0.
#'@param threads Number of threads to do the parallelization. Default is 1.
#'@param savegenenames For which genes their concrete FPM value for each bp 
#'  position need to be saved, or plotted.
#'@param plotgenenames Whether to plot the FPM value for each bp position for 
#'  the genes provided by the parameter \code{savegenenames}. 
#'@param plotfig Whether plot the metagenes or not. Default value is TRUE. For 
#'  the genes indicated by \code{savegenenames}, only if both the paramters 
#'  \code{plotgenenames} and \code{plotfig} are TRUE, they will be plotted.
#'@param genelencutoff The cutoff on gene length (bp). The default value is 
#'  NULL, but if it is set, only genes longer than this cutoff will be 
#'  considerred.
#'@param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM 
#'  value greater than the cutoff will be considerred. Default is 1. 
#'@return Will generate several metagene plots as well as a list containing 
#'  the information of the FPKM values on specific regions for each gene, the 
#'  concrete FPM value on each bp position for the metagene, and the genes 
#'  indicated by the parameter \code{savegenenames}, etc. 
#'@export
metaplot <- function(metafile,
                     targetgenefile, genomename = NULL,
                     tssradius = c(2000, 1000, 500), ttsradius = c(2000, 1000, 500),
                     genebodylen = 4000,
                     label = NULL,
                     strandmethod = 1,
                     threads = 1,
                     savegenenames = NULL, plotgenenames = TRUE, plotfig = TRUE,
                     genelencutoff = NULL,
                     fpkmcutoff = 1){
  
  #Libraries required
  #library(GenomicAlignments)
  
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
      
      genes <- GenomicRanges::GRanges(seqnames = genes$seqnames,
                                      ranges = IRanges::IRanges(start = genes$start, end = genes$end),
                                      strand = genes$strand,
                                      gene_id = genes$gene_id)
      genes <- GenomicAlignments::sort(genes)
      
    }
    
    
  }else if(!is.null(genomename)){
    
    genes <- get(paste0(genomename, '.packagegenes'))
    genes <- GenomicRanges::GRanges(genes)
    
    genes <- genes[genes@ranges@width > max(c(tssradius, ttsradius, genelencutoff))]
    
    genes <- genes[genes$updis >= max(c(tssradius))]
    genes <- genes[genes$dndis >= max(c(ttsradius))]
    
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
  
  #Read in the strand-specific paired-end data#####
  
  reads <- parsereadsfile(readsfile = metafile, strandmethod = strandmethod)
  readsdepth <- length(reads) #FRAGMENT counts
  
  #FPKM screen###
  if(!is.null(fpkmcutoff)){
    
    genes <- getfpkm(features = genes, reads = reads, readsdepth = readsdepth)
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
    
    starts <- subgeneranges@ranges@start - max(tssradius, ttsradius)
    starts[starts < 1] <- 1
    
    widthes <- subgeneranges@ranges@width + max(tssradius, ttsradius)
    ends <- starts + widthes - 1
    endlimits <- reads@seqinfo@seqlengths
    names(endlimits) <- reads@seqinfo@seqnames
    endlimits <- endlimits[as.character(subgeneranges@seqnames)]
    ends[ends > endlimits] <- endlimits[ends > endlimits]
    widthes <- ends - starts + 1
    
    subgeneranges <- GenomicRanges::GRanges(seqnames = subgeneranges@seqnames,
                                            ranges = IRanges::IRanges(start = starts, width = widthes), strand = '*')
    #Set strand as '*' here, so that both the reads on sense strand and antisense strand of a gene
    #can be selected next to plot metagene covering both sense and antisense strands
    
    subreads <- IRanges::subsetByOverlaps(reads, subgeneranges)
    #Here, directly select reads from the paired fragment alignments, no need to separately from first reads
    #alignments and last reads alignments, and convert the strand symbols for one of them from artifical
    #symbols to true symbols, because if select from paired fragment alignments, the strand symbols of these
    #fragments are from the reads alignments file with the true strand symbols, and the one with artifical
    #symbols is disregarded automatically.
    
    subreads <- GenomicAlignments::sort(subreads)
    
    subreadsranges <- range(subreads)
    subgenes <- IRanges::subsetByOverlaps(subgenes, subreadsranges)
    subgenes <- GenomicAlignments::sort(subgenes)
    
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
  
  if(threads == 1){
    
    subreports <- list()
    
    subreports[[1]] <- mapfunction(datalist = parallellist[[1]], depth = readsdepth, 
                                   tssradius = tssradius, ttsradius = ttsradius, 
                                   genebodylen = genebodylen, 
                                   saveindividual = savegenenames)
    
    
  }else{
    
    #library(doParallel)
    #threads <- 8
    
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(min(threads, cores))
    doParallel::registerDoParallel(cl)
    
    date()
    `%dopar%` <- foreach::`%dopar%`
    subreports <- foreach::foreach(datalist = parallellist,
                                   .export = ls(name = globalenv())) %dopar% mapfunction(datalist,
                                                                                         depth = readsdepth,
                                                                                         tssradius = tssradius, ttsradius = ttsradius,
                                                                                         genebodylen = genebodylen,
                                                                                         saveindividual = savegenenames)
    date()
    
    parallel::stopCluster(cl)
    
  }
  
  
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
        
        combinesinglefwdFPMs[[plotsinglegenename]] <- S4Vectors::Rle(singlegenefwdcombine)
        combinesinglerevFPMs[[plotsinglegenename]] <- S4Vectors::Rle(singlegenerevcombine)
        
        
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



#'Generate metagene plots for multiple bam files
#'
#'Generate metagene plots for multiple bam files, as well as other information 
#'including the FPM value of each bp position in the metagene plots, FPKM 
#'value in specific TSS, TTS, or gene body regions of each gene, etc. 
#'
#'@param metafiles The bam files need to generate metagene plots. Should be a 
#'  vector with elements as strings indicating the directores of the files.
#'@param targetgenefile The genes whose FPM values need to be merged together 
#'  to generate metagene plots. If it is NULL, all the genes in the genome 
#'  specified by the parameter \code{genomename} will be analyzed. If provied 
#'  by the user, columns named as chr, start, end, strand, and gene_id are 
#'  required. All genes should be longer than the maximum value of the 
#'  parameters \code{tssradius}, \code{ttsradius}, and \code{genelencutoff}.
#'@param genomename Specify the genome of the genes to be analyzed, when the 
#'  parameter \code{targetgenefile} is NULL.
#'@param tssradius If want to plot the metagene in TSS region, set this 
#'  parameter as a vector with each element corresponding to a radius of the 
#'  TSS region centering around the TSS site. Then, for each radius value, a 
#'  metagene plot will be generated.
#'@param ttsradius If want to plot the metagene in TTS region, set this 
#'  parameter as a vector with each element corresponding to a radius of the 
#'  TTS region centering around the TTS site. Then, for each radius value, a 
#'  metagene plot will be generated. 
#'@param genebodylen If want to plot the metagene in gene body region, set 
#'  this parameter as a specific numeric value (bp). Because different genes 
#'  have different gene body length, before merge their FPM values to a 
#'  metagene, their gene body length will be scaled to a uniform length, which 
#'  is set by this parameter \code{genebodylen}. Genes with a longer length 
#'  will be compressed to this gene length, while the ones with a shorter 
#'  length will be exteneded. 
#'@param labels A vector with elements as strings to be included in the titles 
#'  of the metagene plots to indicate the experimental conditions of the bam 
#'  files. There should be no replicated elements in this vector.
#'@param strandmethod The strand specific method used when preparing the 
#'  sequencing library, can be 1 for directional ligation method and 2 for 
#'  dUTP method. If the sample is sequenced using a single strand method, set 
#'  it as 0.
#'@param threads Number of threads to do the parallelization. Default is 1.
#'@param savegenenames For which genes their concrete FPM value for each bp 
#'  position need to be saved, or plotted.
#'@param plotgenenames Whether to plot the FPM value for each bp position for 
#'  the genes provided by the parameter \code{savegenenames}.
#'@param mergecases Whether merge the data of all bam files together to one, 
#'  and then use it to generate the metagene plots. Default is FALSE.
#'@param genelencutoff The cutoff on gene length (bp). The default value is 
#'  NULL, but if it is set, only genes longer than this cutoff will be 
#'  considerred.
#'@param fpkmcutoff The cutoff on gene FPKM value. Only genes with an FPKM 
#'  value greater than the cutoff will be considerred. Default is 1. 
#'@return Will generate several metagene plots as well as a list with several 
#'  sub-lists containing the information of the FPKM values on specific 
#'  regions for each gene, the concrete FPM value on each bp position for the 
#'  metagene, and the genes indicated by the parameter \code{savegenenames}, 
#'  etc. 
#'@export
mmetaplot <- function(metafiles,
                      targetgenefile,
                      genomename = NULL,
                      tssradius = c(2000, 1000, 500), ttsradius = c(2000, 1000, 500),
                      genebodylen = 4000,
                      labels,
                      strandmethod = 1,
                      threads = 1,
                      savegenenames = NULL, plotgenenames = TRUE,
                      mergecases = FALSE,
                      genelencutoff = NULL,
                      fpkmcutoff = 1){
  
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
      
      mTSSFPKMstatses <- merge(mTSSFPKMstatses, TSSFPKMstatses, 
                               by = c('gene_id'))
      mTTSFPKMstatses <- merge(mTTSFPKMstatses, TTSFPKMstatses, 
                               by = c('gene_id'))
      mGeneBodyFPKMstatses <- merge(mGeneBodyFPKMstatses, GeneBodyFPKMstatses, 
                                    by = c('gene_id'))
      
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
                endtype = c('TSS', 'TTS'), label = NULL)
    
    
    
    mcombinefwdFPMmeans <- as.data.frame(mcombinefwdFPMmeans, stringsAsFactors = FALSE)
    mcombinerevFPMmeans <- as.data.frame(mcombinerevFPMmeans, stringsAsFactors = FALSE)
    colnames(mcombinefwdFPMmeans) <- colnames(mcombinerevFPMmeans) <- names(reslist)
    rownames(mcombinefwdFPMmeans) <- rownames(mcombinerevFPMmeans) <- 1:nrow(mcombinefwdFPMmeans)
    
    plotregions(fwdreads = mcombinefwdFPMmeans, revreads = mcombinerevFPMmeans,
                upwidthlimit = max(tssradius), upwidth = NULL,
                dnwidthlimit = max(ttsradius), dnwidth = NULL,
                endtype = c('TSS', 'TTS'), label = NULL)
    
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


