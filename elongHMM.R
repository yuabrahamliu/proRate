######
##
##Created on Jun 5, 2018
##
##@author: Yu Liu
##
######

library(GenomicAlignments)
library(depmixS4)

parseConfig <- function(){
  configs <- readLines('config')
  idx <- grep('time2file =', configs)
  time2file <- configs[idx]
  idx <- grep('samplepair =', configs)
  samplepair <- configs[idx]
  idx <- grep('time1file =', configs)
  time1file <- configs[idx]
  idx <- grep('targetfile =', configs)
  targetfile <- configs[idx]
  idx <- grep('time =', configs)
  time <- configs[idx]
  idx <- grep('strandmethod =', configs)
  strandmethod <- configs[idx]
  idx <- grep('kscutoff =', configs)
  kscutoff <- configs[idx]
  idx <- grep('prirate = ', configs)
  prirate <- configs[idx]
  
  time2file <- unlist(strsplit(time2file, split = ' = '))[2]
  if(is.na(time2file)){
    cat('Please provide the time point2 bam file\n')
    q(save = 'no')
  }
  samplepair <- unlist(strsplit(samplepair, split = ' = '))[2]
  if(is.na(samplepair)){
    cat('Please provide the sample pair name\n')
    q(save = 'no')
  }
  time1file <- unlist(strsplit(time1file, split = ' = '))[2]
  if(is.na(time1file)){
    cat('Please provide the time point1 bam file\n')
    q(save = 'no')
  }
  targetfile <- unlist(strsplit(targetfile, split = ' = '))[2]
  if(is.na(targetfile)){
    cat('Please provide the target gene list\n')
    q(save = 'no')
  }
  time <- as.numeric(unlist(strsplit(time, split = ' = '))[2])
  if(is.na(time)){
    cat('Please provide the elongation time\n')
    q(save = 'no')
  }
  strandmethod <- unlist(strsplit(strandmethod, split = ' = '))[2]
  if(is.na(strandmethod)){
    strandmethod <- 1
  }else if(strandmethod != 'single'){
    strandmethod <- as.numeric(strandmethod)
  }
  kscutoff <- as.numeric(unlist(strsplit(kscutoff, split = ' = '))[2])
  if(is.na(kscutoff)){
    kscutoff <- 0.01
  }
  prirate <- as.numeric(unlist(strsplit(prirate, split = ' = '))[2])
  if(is.na(prirate)){
    prirate <- 'None'
  }
  return(list(time2file, samplepair, time1file, targetfile, 
              time, strandmethod, kscutoff, prirate))

}

windowAnalysis <- function(reads, strand="*", chrom=NULL) {
  
  if (!is.null(chrom))
    reads <- reads[seqnames(reads) == chrom,]
  
  readsList <- split(reads, seqnames(reads))
  #important!
  #'split' divides the data in a vector-like object 'x' into the
  #groups defined by 'f'.
  #Here split the bam file data into different chroms
  
  # Change reads strand 
  readsList <- endoapply(readsList, function(x) {
    if (strand == "*") 
      strand(x) <- "*"
    else
      x <- x[strand(x) == strand,]
    x 
  })
  
  #extract the specific strand reads defined by the fuction parameter
  
  #lapply returns a list of the same length as X, each element of which is 
  #the result of applying FUN to the corresponding element of X.
  
  #sapply is a user-friendly version and wrapper of lapply by default 
  #returning a vector, matrix.
  
  #endoapply performs the endomorphic equivalents of lapply 
  #by returning objects of the same class as the inputs rather than a list.
  
  testlist <- as.list(readsList)
  
  if(!is.null(chrom)){
    testlist <- testlist[chrom]
  }else{
    chroms <- names(testlist)
    normchroms <- chroms[grepl('chr', chroms)]
    normchroms <- normchroms[normchroms != 'chrM']
    normchroms <- unique(normchroms)
    testlist <- testlist[normchroms]
    normchroms <- names(testlist)[lapply(testlist, length) > 0]
    testlist <- testlist[normchroms]
  }
  
  H <- list()
  for(i in 1:length(testlist)){
    x <- testlist[[i]]
    seqlevels(x) <- seqlevelsInUse(x)
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

rateHMM <- function(){
  ###################
  #Config
  params <- parseConfig()
  time2file <- params[[1]]
  samplepair <- params[[2]]
  time1file <- params[[3]]
  targetfile <- params[[4]]
  txtname <- paste0(samplepair, '_wave.txt')
  time <- params[[5]]
  strandmethod <- params[[6]]
  kscutoff <- params[[7]]
  prirate <- params[[8]]
  
  ###################
  setwd(paste0(getwd(), '/elongHMM_pipeline/data'))
  #Prepare the genes to GRanges object
  genes <- read.table(targetfile, sep = '\t', header = TRUE, 
                      stringsAsFactors = FALSE, quote = '', 
                      check.names = FALSE)
  genes$start <- genes$start + 1
  names(genes) <- c('seqnames', 'start', 'end', 'strand', 'SYMBOL', 'ID')
  
  ###################
  #Read in the strand-specific paired-end data
  if(strandmethod != 'single'){
    reads1p_ori <- readGAlignmentPairs(time1file, strandMode = strandmethod)
  }else{
    reads1p_ori <- readGAlignments(time1file)
  }
  reads1 <- as(reads1p_ori, 'GRanges')
  rm(reads1p_ori)
  
  if(strandmethod != 'single'){
    reads2p_ori <- readGAlignmentPairs(time2file, strandMode = strandmethod)
  }else{
    reads2p_ori <- readGAlignments(time2file)
  }
  reads2 <- as(reads2p_ori, 'GRanges')
  rm(reads2p_ori)
  ###################
  #Other default parameters
  window_num <- 40
  ###################
  #Convert the 2 bams file to strand separated wig files 
  Fp1 <- windowAnalysis(reads=reads1, strand="+")
  Fp2 <- windowAnalysis(reads=reads2, strand="+")
  Fm1 <- windowAnalysis(reads=reads1, strand="-")
  Fm2 <- windowAnalysis(reads=reads2, strand="-")
  ###################
  #Number of states in HMM!
  nstates<- 2
  startshorten <- 1000
  endshorten <- 1000
  
  syms <- c()
  distances <- c()
  times <- c()
  rates <- c()
  
  for(i in 1:nrow(genes)){
    ###################
    print(i)
    gene_len <- abs(genes[i, 3] - genes[i, 2]) + 1
    
    if(gene_len < 2*(startshorten + endshorten)){
      next()
    }
    
    if(genes[i,4] == "+") {
      
      start <- genes[i, 2] + startshorten
      end   <- genes[i, 3] - endshorten
      
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
        (as.numeric(Fp1[[ as.character(genes[i,1])]]))[c(start:end)]
      emis2  <- 
        (as.numeric(Fp2[[ as.character(genes[i,1])]]))[c(start:end)]
    }else {
      start <- genes[i, 2] + endshorten
      end   <- genes[i, 3] - startshorten
      emis1  <- 
        rev((as.integer(Fm1[[ as.character(genes[i,1])]]))
            [c(start:end)])
      #Note here the minus strand has been reversed!
      emis2  <- 
        rev((as.integer(Fm2[[ as.character(genes[i,1])]]))
            [c(start:end)])
      #Note here the minus strand has been reversed!
    }
    
    depth1 <- sum(emis1)
    depth2 <- sum(emis2)
    
    if(depth1 == 0 & depth2 == 0){
      next()
    }
    
    if(depth1 == 0 | depth2 == 0){
      emis1 <- emis1 + 10^-5
      emis2 <- emis2 + 10^-5
      depth1 <- sum(emis1)
      depth2 <- sum(emis2)
    }
    
    weight1 <- depth1/mean(c(depth1, depth2))
    weight2 <- depth2/mean(c(depth1, depth2))
    
    emis1_adj <- emis1/weight1
    emis2_adj <- emis2/weight2
    
    #Use ks.test to primarily eliminate genes without significant wave change
    ksresult <- ks.test(emis1_adj, emis2_adj)
    kspval <- ksresult$p.val
    if(kspval >= kscutoff){
      next()
    }
    
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
    if(sum(is.na(gene)) > 0 | sum(is.infinite(gene)) > 0){
      windows1 <- windows1 + 10^-5
      windows2 <- windows2 + 10^-5
      gene <- windows2/windows1
    }
    
    ###################
    #Begin to infer
    gene <- data.frame(ratio_adj = gene, stringsAsFactors = FALSE)
    
    set.seed(1)
    if(prirate == 'None'){
      mod <- depmix(response = ratio_adj ~ 1, data = gene, nstates = nstates, 
                    instart = c(1, 0))
    }else{
      approxDist <- prirate*time
      mod <- depmix(response = ratio_adj ~ 1, data = gene, nstates = nstates, 
                    instart = c(1, 0), 
                    trstart = c(round((approxDist-startshorten)/size)/window_num, 
                                (1 - round((approxDist-startshorten)/size)/window_num), 
                                0.00001, 0.99999))
      rm(approxDist)
    }
    
    fm <- tryCatch({
      fit(mod)
    }, error = function(err){
      NULL
    })
    
    if(is.null(fm)){
      next()
    }else{
      fmpost <- posterior(fm)
    }
    
    states <- fmpost$state
    #select transition point
    if(length(unique(states)) == 1){
      next()
    }
    
    meandiff_state <- c()
    for(j in 1:(length(states)-1)) {
      left_state  <- states[c(1:j)]
      right_state <- states[c((j+1):length(states))]
      left_state_mean <- mean(left_state)
      right_state_mean <- mean(right_state)
      diff_state_rl <- right_state_mean - left_state_mean
      meandiff_state <- c(meandiff_state, diff_state_rl)
    }
    
    point_idx <- match(max(meandiff_state), meandiff_state)
    
    if(is.na(point_idx)){
      next()
    }

    ###################
    forward <- round(size * (point_idx - 1))
    backward <- round(size * (point_idx + 1))
    
    #focus on the single window with transition point
    if(gene[point_idx, 1] == gene[(point_idx + 1), 1]){
      final_point <- size*point_idx
    }else{
      
      startcoord <- 1 + forward
      endcoord <- backward
      
      expand_ratio <- (emis2_adj[startcoord:endcoord])/
        (emis1_adj[startcoord:endcoord])
      if(sum(is.na(expand_ratio)) > 0 | sum(is.infinite(expand_ratio)) > 0){
        expand_ratio <- (emis2_adj[startcoord:endcoord] + 10^-5)/
          (emis1_adj[startcoord:endcoord] + 10^-5)
      }
      
      square_list <- c()
      for(k in startcoord:(endcoord - 1)){
        leftidx <- startcoord:k
        rightidx <- (k+1):endcoord
        leftratio <- expand_ratio[leftidx - startcoord + 1]
        rightratio <- expand_ratio[rightidx - startcoord + 1]
        
        leftratio_var <- sum((leftratio - mean(leftratio))^2)
        rightratio_var <- sum((rightratio - mean(rightratio))^2)
        square <- leftratio_var + rightratio_var
        square_list <- c(square_list, square)
      }
      
      final_point <- match(min(square_list), square_list)
    }
    
    distance <- startshorten + forward + final_point
    
    syms <- c(syms, genes[i, 5])
    distances <- c(distances, distance)
    times <- c(times, time)
    rates <- c(rates, distance/time)
  }
  
  report <- data.frame(genesym = syms, 
                       distance = distances, 
                       time = times, 
                       rate = rates, 
                       stringsAsFactors = FALSE)
  setwd('../result')
  write.table(report, txtname, sep = '\t', 
              quote = FALSE, row.names = FALSE)

}

#####

if(!interactive()){
  rateHMM()
}













