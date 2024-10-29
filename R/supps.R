
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
    vi <- IRanges::Views(cov, start=starts, width=1)
    
    H[[i]] <- S4Vectors::Rle(IRanges::viewSums(vi))
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



#'Process the results of the functions \code{metaplot} and \code{mmetaplot} 
#'
#'Process the metagene results from \code{metaplot} and \code{mmetaplot} and 
#'remove the extremely large FPM outliers in the original metagene plots.
#'
#'@param fwdlist A list containing the metagene FPM values from the function 
#'  \code{metaplot} or \code{mmetaplot}'s result list. Each element in it is 
#'  an FPM value vector contained in the original result list, and its slot 
#'  name there is "TSSfwdFPMmeans", "TTSfwdFPMmeans", "GeneBodyfwdFPMmeans", 
#'  or "combinefwdFPMmeans". It is the FPM value vector for the metagene's 
#'  sense strand and should match the corresponding list element transferred 
#'  to another parameter \code{revlist}, which is the FPM value vector for the 
#'  metagene's antisense strand. For the different elements of \code{fwdlist}, 
#'  they should come from the same regions, such as the TSS region, the TTS 
#'  region, the gene body region, etc. Their difference is that they may from 
#'  different experimental conditions, such as the TSS region's metagene data 
#'  of wildtype samples and gene knockout samples, respectively.
#'@param revlist A list containing the metagene FPM values from the function 
#'  \code{metaplot} or \code{mmetaplot}'s result list. Each element in it is 
#'  an FPM value vector contained in the original result list, and its slot 
#'  name there is "TSSrevFPMmeans", "TTSrevFPMmeans", "GeneBodyrevFPMmeans", 
#'  or "combinerevFPMmeans". It is the FPM value vector for the metagene's 
#'  antisense strand and should match the list element transferred to another 
#'  parameter \code{fwdlist}, which is the FPM value vector for the metagene's 
#'  sense strand. For the different elements of \code{fwdlist}, they should 
#'  come from the same regions, such as the TSS region, the TTS region, the 
#'  gene body region, etc. Their difference is that they may from different 
#'  experimental conditions, such as the TSS region's metagene antisense FPM 
#'  values of wildtype samples and gene knockout samples, respectively.
#'@param groupnames A vector with elements as strings to show the different 
#'  elements' experimental conditions in \code{fwdlist} and \code{revlis}.
#'@param lineposes A vector with elements as integers to show the coordinates 
#'  of the TSS and/or TTS points in the metagene. These coordinates use the 
#'  beginning of the metagene x-coordinate as 1.
#'@param labels A vector with elements as strings to show the x-axis labels 
#'  for the beginning and end of the metagene's x-axis and for the positions 
#'  indicated by the parameter \code{lineposes}. The labels should be ordered 
#'  following their coordinates along the x-axis.
#'@param title The title for the metagene plot. 
#'@param cutoff To remove the extremely large FPM outliers from the original 
#'  metagene plot, this parameter needs to be set. The default value is 0.01, 
#'  which means the 1 - 1% (0.01) = 99% quantile of the metagene FPMs will be 
#'  defined as their maximum, and any larger ones will be reduced. 
#'@param titlesize The font size for the plot title. Default is 17. 
#'@param textsize The font size for the plot texts. Default is 16.
#'@return The processed metagene plot with the FPM outliers removed. 
#'
#'
#'
#'@examples
#'library(proRate)
#'
#'wt0file <- system.file("extdata", "wt0.bam", package = "proRate")
#'ko0file <- system.file("extdata", "ko0.bam", package = "proRate")
#'
#'metareslist <- mmetaplot(metafiles = c(wt0file, ko0file), 
#'                         labels = c("WT", "KO"), 
#'                         tssradius = c(1000, 500), 
#'                         ttsradius = c(1000), 
#'                         genebodylen = 2000, 
#'                         strandmethod = 1, 
#'                         genomename = "mm10", 
#'                         genelencutoff = 40000, 
#'                         fpkmcutoff = 1)
#'
#'combinefwdlist <- list()
#'
#'combinerevlist <- list()
#'
#'for(i in seq(1, 2, 1)){
#'  
#'  groupname <- c("WT", "KO")[i]
#'  
#'  combinefwdlist[[i]] <- metareslist[[groupname]]$combinefwdFPMmeans
#'  
#'  combinerevlist[[i]] <- metareslist[[groupname]]$combinerevFPMmeans
#'  
#'}
#'
#'plotprocessing(fwdlist = combinefwdlist, 
#'               revlist = combinerevlist, 
#'               
#'               cutoff = 0.01, 
#'               groupnames = c("WT", "KO"), 
#'               labels = c("-1000", "TSS", "TTS", "+1000"), 
#'               lineposes = c(1001, 3000), 
#'               
#'               title = "WT_KO metagene from -1000bp of TSS to +1000bp of TTS", 
#'               titlesize = 17, 
#'               textsize = 16)
#'
#'
#'
#'@export
plotprocessing <- function(fwdlist, 
                           revlist, 
                           groupnames, 
                           lineposes, 
                           labels, 
                           title, 
                           cutoff = 0.01, 
                           titlesize = 17, 
                           textsize = 16){
  
  reads_totals <- list()
  
  i <- 1
  
  for(i in 1:length(groupnames)){
    
    groupname <- groupnames[i]
    
    reads_total <- data.frame(xcoord = c(seq(1, length(fwdlist[[i]]), 1), 
                                         seq(1, length(revlist[[i]]), 1)), 
                              reads = c(fwdlist[[i]], -revlist[[i]]), 
                              condition = groupname, 
                              strand = rep(c('Sense', 'Antisense'), 
                                           c(length(fwdlist[[i]]), 
                                             length(revlist[[i]]))), 
                              stringsAsFactors = FALSE)
    
    reads_totals[[i]] <- reads_total
    
  }
  
  reads_totals <- do.call(rbind, reads_totals)
  
  reads_totals$condition <- factor(reads_totals$condition, 
                                   levels = groupnames, 
                                   ordered = TRUE)
  
  reads_totals$strand <- factor(reads_totals$strand, 
                                levels = c('Sense', 'Antisense'), 
                                ordered = TRUE)
  
  if(!is.null(cutoff)){
    
    sensecutoff <- quantile(subset(reads_totals, strand == 'Sense')$reads, 
                            1 - cutoff)
    
    antisensecutoff <- -quantile(abs(subset(reads_totals, strand == 'Antisense')$reads), 
                                 1 - cutoff)
    
    reads_totals$reads[reads_totals$reads > sensecutoff & 
                         reads_totals$strand == 'Sense'] <- sensecutoff
    
    reads_totals$reads[reads_totals$reads < antisensecutoff & 
                         reads_totals$strand == 'Antisense'] <- antisensecutoff
    
  }
  
  points <- unique(c(1, lineposes, max(reads_totals$xcoord)))
  
  pointlabels <- labels
  
  p <- ggplot2::ggplot(reads_totals, 
                       mapping = ggplot2::aes(x = xcoord, y = reads, color = condition))
  
  p <- p + ggplot2::geom_line(size = 1) + 
    
    ggplot2::scale_x_continuous(breaks = points,
                                labels = pointlabels) + 
    
    ggplot2::xlab('') + ggplot2::ylab('FPM') + 
    
    ggplot2::geom_vline(xintercept = lineposes[1], linetype = 2, color = 'red', size = 1) +
    ggplot2::geom_vline(xintercept = lineposes[length(lineposes)], linetype = 2, color = 'red', size = 1) +
    ggplot2::facet_grid(strand~., scales = 'free_y') +
    ggplot2::ggtitle(title) +
    ggplot2::scale_color_discrete(name = 'Condition') +
    ggplot2::theme_bw() + 
    
    ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = 'bold'),
                   #plot.subtitle = ggplot2::element_text(size = textsize, face = 'bold'),
                   axis.title = ggplot2::element_text(size = titlesize, face = 'bold'),
                   axis.text = ggplot2::element_text(size = textsize, face = 'bold'), 
                   axis.text.x = ggplot2::element_text(hjust=1, angle = 90), 
                   legend.title = ggplot2::element_text(size = textsize, face = 'bold'),
                   legend.text = ggplot2::element_text(size = textsize, face = 'bold'),
                   strip.text = ggplot2::element_text(size = textsize, face = 'bold'), 
                   panel.grid = ggplot2::element_blank(), 
                   panel.spacing = grid::unit(0, 'lines'))
  
  if(length(unique(reads_totals$condition)) == 1){
    p <- p + ggplot2::guides(color = FALSE)
    
  }
  
  print(p)
  
}



#'Add gene coordinate track to a gene's rate inference result plot 
#'
#'Add gene coordinate track to a gene's rate inference result plot generated 
#'by the function \code{calrate} or \code{mcalrate}. 
#'
#'@param genedat The gene rate inference report data frame generated by the 
#'  function \code{calrate} or \code{mcalrate}. It can be extracted from the 
#'  "report" slot of these functions' result, which is a data frame, and what 
#'  needed here is the its sub-data-frame (with only one row) containing the 
#'  result information of the specific gene whose inference plots need to be 
#'  modified here.
#'@param binplotdat The bin-level rate inference plot data of a specific gene. 
#'  It can be extracted from the "binplots" slot of the results generated by 
#'  \code{calrate} or \code{mcalrate}. What needed here is the element data 
#'  contained in this slot using the gene's name as the element name.
#'@param expandplotdat The base-pair-level rate inference plot data of a gene. 
#'  It can be extracted from the "expandplots" slot of the results generated 
#'  by the function \code{calrate} or \code{mcalrate}. What needed here is the 
#'  element data contained in this slot using the gene's name as the element 
#'  name.
#'@param genomename Specify the genome of the specific gene whose inference 
#'  plots need to be modified here. It can be "mm10" for mouse or "hg38" for 
#'  human.
#'@param textsize The font size for the plot texts. Default is 13.
#'@param titlesize The font size for the plot titles. Default is 15.
#'@param face The font face for the plot texts. Default is "bold".
#'@param method The method used when inferring the gene's transcription rate 
#'  with the function \code{calrate} or \code{mcalrate}. Can be "LSS" for the 
#'  least sum of squares method, or "HMM" for the hidden Markov model method.
#'@return The modified gene rate inference plots with the gene's coordinate 
#'  track added.
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
#'                   threads = 4, 
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
#'addtrack(genedat = subset(wtrates$report, gene_id == "Mamdc2"), 
#'         binplotdat = wtrates$binplots$Mamdc2, 
#'         expandplotdat = wtrates$expandplots$Mamdc2, 
#'         genomename = "mm10", 
#'         method = "LSS", 
#'         titlesize = 17, 
#'         textsize = 16, 
#'         face = "bold")
#'
#'
#'
#'@export
addtrack <- function(genedat, 
                     binplotdat, 
                     expandplotdat, 
                     genomename = 'mm10', 
                     textsize = 13, 
                     titlesize = 15, 
                     face = 'bold', 
                     method = NULL){
  
  binplottitle <- binplotdat$labels$title
  
  expandplottitle <- expandplotdat$labels$title
  
  if(!is.null(method)){
    
    binplottitle <- gsub(pattern = '\\(Distance', 
                         replacement = paste0('(', method, ' Distance'), 
                         binplottitle)
    
    expandplottitle <- gsub(pattern = '\\(Distance', 
                            replacement = paste0('(', method, ' Distance'), 
                            expandplottitle)
    
  }
  
  
  
  windownum <- gsub(pattern = '^.*Length = ', replacement = '', x = binplottitle)
  windownum <- gsub(pattern = '\\,.*$', replacement = '', x = windownum)
  windownum <- as.numeric(windownum)
  
  windowdistance <- gsub(pattern = '^.*Distance = ', replacement = '', x = binplottitle)
  windowdistance <- gsub(pattern = '\\,.*$', replacement = '', x = windowdistance)
  windowdistance <- as.numeric(windowdistance)
  
  windowsize <- gsub(pattern = '^.*Length = ', replacement = '', x = expandplottitle)
  windowsize <- gsub(pattern = '\\,.*$', replacement = '', x = windowsize)
  windowsize <- as.numeric(windowsize)/2
  
  expanddistance <- gsub(pattern = '^.*Distance = ', replacement = '', x = expandplottitle)
  expanddistance <- gsub(pattern = '\\,.*$', replacement = '', x = expanddistance)
  expanddistance <- as.numeric(expanddistance)
  
  startshorten <- genedat$distance - windowsize * (windowdistance - 1) - expanddistance
  
  endshorten <- genedat$genewidth - windowsize * windownum - startshorten
  
  expandstart <- startshorten + windowsize * (windowdistance - 1)
  
  expandend <- startshorten + windowsize * (windowdistance + 1) - 1
  
  
  
  points <- c(1, binplotdat$plot_env$distance, binplotdat$plot_env$end)
  pointlabels <- c(binplotdat$plot_env$endtype[1], 'Transition', 
                   binplotdat$plot_env$endtype[2])
  
  expandpoints <- c(1, expanddistance, expandplotdat$plot_env$end)
  expandpointlabels <- c(sub(pattern = '\\+[0-9].*bp$', 
                             replacement = paste0('+', expandstart, 'bp'), 
                             x = expandplotdat$plot_env$endtype[1]), 
                         'Transition', 
                         sub(pattern = '\\+[0-9].*bp$', 
                             replacement = paste0('+', expandend, 'bp'), 
                             x = expandplotdat$plot_env$endtype[2]))
  
  ylabel <- binplotdat$plot_env$ylabel
  
  expandylabel <- expandplotdat$plot_env$ylabel
  
  
  
  chr <- genedat$chr
  start <- genedat$start
  end <- genedat$end
  strand <- genedat$strand
  
  genegr <- GenomicRanges::GRanges(seqnames = chr, 
                                   ranges = IRanges::IRanges(start = start, 
                                                             end = end), 
                                   strand = strand)
  
  genegr$gene_id <- genedat$gene_id
  
  
  
  binplotstart <- start + startshorten
  binplotend <- end - endshorten
  
  binplotgr <- GenomicRanges::GRanges(seqnames = chr, 
                                      ranges = IRanges::IRanges(start = binplotstart, 
                                                                end = binplotend), 
                                      strand = strand)
  
  if(strand != '-'){
    
    expandplotstart <- start + expandstart
    expandplotend <- start + expandend
    
  }else{
    
    expandplotstart <- end - expandend
    expandplotend <- end - expandstart
    
  }
  
  
  
  expandplotgr <- GenomicRanges::GRanges(seqnames = chr, 
                                         ranges = IRanges::IRanges(start = expandplotstart, 
                                                                   end = expandplotend), 
                                         strand = strand)
  
  
  
  if(genomename == 'hg38'){
    
    #library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    #library(org.Hs.eg.db)
    
    exoncoords <- 
      GenomicFeatures::exons(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene)
    
  }else if(genomename == 'mm10'){
    
    #library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    #library(org.Mm.eg.db)
    
    exoncoords <- 
      GenomicFeatures::exons(TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene)
    
  }
  
  geneexons <- GenomicRanges::intersect(exoncoords, genegr)
  
  if(length(geneexons) == 0){
    
    return(NULL)
    
  }else{
    
    geneexons <- GenomicRanges::reduce(geneexons)
    geneexons$gene_id <- genegr$gene_id
    
  }
  
  split_body <- function(df, width = 1000){
    
    wd <- df$end - df$start + 1
    
    strand <- unique(df$strand)
    
    nbreak <- wd/width
    
    if(df$strand == '-'){
      
      if(nbreak > 1){
        
        steps <- 0:nbreak
        starts <- width*steps + df$start
        
        starts[starts > df$end] <- NA
        starts <- starts[!is.na(starts)]
        
      }else{
        
        starts <- 0
        
      }
      
    }else{
      
      if(nbreak > 1){
        
        steps <- 0:nbreak
        starts <- width*steps + df$start
        
        starts[starts >= df$end] <- NA
        starts[steps == 0] <- NA
        starts <- starts[!is.na(starts)]
        starts <- c(starts, df$end)
        
      }else{
        
        starts <- df$end
        
      }
      
    }
    
    breaks <- data.frame(
      seqnames = df$seqnames[[1]],
      start = starts,
      end = starts,
      strand = df$strand[[1]],
      gene_name = df$gene_name[[1]], 
      type = 'arrow', 
      gene_biotype = df$gene_biotype[[1]]
    )
    
    return(breaks)
    
  }
  
  reformat_annotations <- function(annotation, 
                                   start.pos, 
                                   end.pos){
    
    total.width <- end.pos - start.pos + 1
    tick.freq <- total.width / 50
    exons <- as.data.frame(annotation)
    
    annotation <- as.data.frame(annotation)
    
    # add gene total start / end
    gene_bodies <- list()
    
    
    
    df <- data.frame(seqnames = annotation$seqnames[[1]], 
                     start = min(annotation$start), 
                     end = max(annotation$end), 
                     strand = annotation$strand[[1]], 
                     gene_name = annotation$gene_id[[1]], 
                     type = 'body', 
                     gene_biotype = character(1))
    
    # trim any that extend beyond region
    df$start <- ifelse(
      
      test = df$start < start.pos,
      yes = start.pos,
      no = df$start
      
    )
    
    df$start <- start.pos
    
    df$end <- ifelse(
      
      test = df$end > end.pos,
      yes = end.pos,
      no = df$end
      
    )
    
    df$end <- end.pos
    
    breaks <- split_body(df = df, width = tick.freq)
    
    df <- rbind(df, breaks)
    
    
    
    gene_bodies <- df
    
    label_df <- gene_bodies[gene_bodies$type == 'body',]
    label_df$width <- label_df$end - label_df$start + 1
    label_df$position <- label_df$start + (label_df$width / 2)
    
    onplus <- gene_bodies[gene_bodies$strand %in% c('*', '+') & 
                            gene_bodies$type == 'arrow', ]
    onminus <- gene_bodies[gene_bodies$strand == '-' & 
                             gene_bodies$type == 'arrow', ]
    
    reslist <- list()
    
    reslist$labels <- label_df
    reslist$exons <- exons
    reslist$plus <- onplus
    reslist$minus <- onminus
    
    return(reslist)
    
  }
  
  geneplot <- function(region, 
                       geneexons, 
                       genegr){
    
    annotation.subset <- IRanges::subsetByOverlaps(x = geneexons, ranges = region)
    
    if(length(annotation.subset) > 0){
      
      annotation.subset.start <- min(annotation.subset@ranges@start)
      
      annotation.subset.end <- max(annotation.subset@ranges@start + 
                                     annotation.subset@ranges@width - 1)
      
      geneexons.starts <- geneexons@ranges@start
      
      geneexons.ends <- geneexons@ranges@start + geneexons@ranges@width - 1
      
      geneexons.start <- min(geneexons.starts)
      
      geneexons.end <- max(geneexons.ends)
      
      
      region.start <- region@ranges@start
      
      region.end <- as.integer(region@ranges@start + region@ranges@width - 1)
      
      geneplotstart <- region.start
      
      geneplotend <- region.end
      
      
      
      if(any(annotation.subset.start == geneexons.starts) & 
         (annotation.subset.start < region.start)){
        
        
        
        
        
        annotation.subset@ranges@width[annotation.subset@ranges@start == annotation.subset.start] <- 
          annotation.subset@ranges@width[annotation.subset@ranges@start == annotation.subset.start] - 
          (region.start - annotation.subset.start)
        
        
        
        
        
        annotation.subset@ranges@start[annotation.subset@ranges@start == annotation.subset.start] <- 
          region.start
        
      }
      
      if(any(annotation.subset.end == geneexons.ends) & 
         (annotation.subset.end > region.end)){
        
        annotation.subset@ranges@width[annotation.subset@ranges@start + annotation.subset@ranges@width - 1 == 
                                         annotation.subset.end] <- 
          region.end - annotation.subset@ranges@start[annotation.subset@ranges@start + annotation.subset@ranges@width - 1 == 
                                                        annotation.subset.end] + as.integer(1)
        
      }
      
      if(any(annotation.subset.start == geneexons.starts) & 
         (annotation.subset.start > region.start)){
        
        #geneplotstart <- annotation.subset.start
        
      }
      
      if(any(annotation.subset.end == geneexons.ends) & 
         (annotation.subset.end < region.end)){
        
        #geneplotend <- annotation.subset.end
        
      }
      
      annotation_df_list <- reformat_annotations(
        
        annotation = annotation.subset,
        start.pos = geneplotstart,
        end.pos = geneplotend
        
      )
      
    }else{
      
      annotation_df_list <- list()
      
      annotation_df_list$labels <- data.frame(seqnames = geneexons@seqnames@values, 
                                              start = region@ranges@start, 
                                              end = region@ranges@start + region@ranges@width - 1, 
                                              strand = region@strand@values, 
                                              gene_name = unique(genegr$gene_id), 
                                              type = 'body', 
                                              gene_biotype = character(1), 
                                              width = region@ranges@width, 
                                              position = region@ranges@start + 
                                                (region@ranges@width / 2))
      
      annotation_df_list$exons <- data.frame(seqnames = factor(x = character(0), 
                                                               levels(annotation_df_list$labels$seqnames)), 
                                             start = integer(0), 
                                             end = integer(0), 
                                             width = integer(0), 
                                             strand = factor(x = character(0), 
                                                             levels(annotation_df_list$labels$strand)), 
                                             gene_id = character(0))
      
      df <- data.frame(seqnames = annotation_df_list$labels$seqnames, 
                       start = annotation_df_list$labels$start, 
                       end = annotation_df_list$labels$end, 
                       strand = annotation_df_list$labels$strand, 
                       gene_name = annotation_df_list$labels$gene_name, 
                       type = 'body', 
                       gene_biotype = character(1))
      
      total.width <- annotation_df_list$labels$end - annotation_df_list$labels$start + 1
      
      tick.freq <- total.width / 50
      
      breaks <- split_body(df = df, width = tick.freq)
      
      df <- rbind(df, breaks)
      
      onplus <- df[df$strand %in% c('*', '+') & 
                     df$type == 'arrow', ]
      onminus <- df[df$strand == '-' & 
                      df$type == 'arrow', ]
      
      annotation_df_list$plus <- onplus
      
      annotation_df_list$minus <- onminus
      
    }
    
    annotation_df_list$geneplotlabel <- 'Complete'
    
    if(annotation_df_list$labels$start > genegr@ranges@start){
      
      annotation_df_list$geneplotlabel <- 'Truncated'
      
    }
    
    if(annotation_df_list$labels$end < genegr@ranges@start + genegr@ranges@width - 1){
      
      annotation_df_list$geneplotlabel <- 'Truncated'
      
    }
    
    return(annotation_df_list)
    
  }
  
  annotation_df_list <- geneplot(region = genegr, geneexons = geneexons, 
                                 genegr = genegr)
  
  #annotation_df_list <- geneplot(region = binplotgr, geneexons = geneexons, 
  #                               genegr = genegr)
  
  expandannotation_df_list <- geneplot(region = expandplotgr, geneexons = geneexons, 
                                       genegr = genegr)
  
  
  
  reads_totals <- binplotdat$data
  
  expandreads_totals <- expandplotdat$data
  
  if(genegr@strand@values != '-'){
    
    reads_totals$xcoord <- start + startshorten + 
      ceiling(reads_totals$xcoord * windowsize)
    
    points <- start + startshorten + 
      ceiling(points * windowsize)
    
    expandreads_totals$xcoord <- seq(expandplotstart, expandplotend, 1)
    
    expandpoints <- c(start + expandstart, 
                      start + expandstart + expandpoints[2] - 1, 
                      start + expandend)
    
  }else{
    
    reads_totals$xcoord <- end - startshorten + 1 - 
      floor(reads_totals$xcoord * windowsize)
    
    reads_totals <- reads_totals[order(-seq(row.names(reads_totals))), , drop = TRUE]
    
    row.names(reads_totals) <- 1:nrow(reads_totals)
    
    points <- end - startshorten + 1 - 
      floor(points * windowsize)
    
    expandreads_totals$xcoord <- rev(seq(expandplotstart, expandplotend, 1))
    
    expandreads_totals <- expandreads_totals[order(-seq(row.names(expandreads_totals))), , 
                                             drop = TRUE]
    
    row.names(expandreads_totals) <- 1:nrow(expandreads_totals)
    
    expandpoints <- c(end - expandend, 
                      end - expandstart - expandpoints[2] + 1, 
                      end - expandstart)
    
  }
  
  
  
  p <- ggplot2::ggplot(reads_totals, mapping = ggplot2::aes(x = xcoord, y = reads))
  
  p <- p + ggplot2::geom_point(size = 2, color = scales::hue_pal()(1))
  
  ylims <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$continuous_range
  
  xlabel <- paste0(genedat$gene_id, ' ', chr, 
                   ':', annotation_df_list$labels$start, '-', 
                   annotation_df_list$labels$end, 'bp')
  
  if(annotation_df_list$geneplotlabel == 'Truncated'){
    
    xlabel <- paste0(genedat$gene_id, ' (Truncated) ', chr, 
                     ':', annotation_df_list$labels$start, '-', 
                     annotation_df_list$labels$end, 'bp')
    
  }
  
  p <- p +
    # exons
    ggplot2::geom_segment(
      
      data = annotation_df_list$exons,
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        xend = end,
        yend = ylims[1],
        color = strand
        
      ), 
      
      show.legend = FALSE,
      linewidth = 6
      
    ) + 
    # gene body
    ggplot2::geom_segment(
      
      data = annotation_df_list$labels,
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        xend = end,
        yend = ylims[1],
        color = strand
        
      ), 
      
      show.legend = FALSE,
      linewidth = 1
      
    )
  
  if(nrow(annotation_df_list$plus) > 0){
    # forward strand arrows
    p <- p + ggplot2::geom_segment(
      
      data = annotation_df_list$plus,
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        xend = end,
        yend = ylims[1],
        color = strand
        
      ), 
      
      arrow = ggplot2::arrow(
        
        ends = 'first',
        type = 'open',
        angle = 135,
        length = ggplot2::unit(x = 0.075, units = 'inches')
        
      ), 
      
      show.legend = FALSE,
      linewidth = 1/2
      
    )
    
  }
  
  if(nrow(annotation_df_list$minus) > 0){
    
    # reverse strand arrows
    p <- p + ggplot2::geom_segment(
      
      data = annotation_df_list$minus, 
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        xend = end,
        yend = ylims[1],
        color = strand
        
      ), 
      
      arrow = ggplot2::arrow(
        
        ends = 'first',
        type = 'open',
        angle = 45,
        length = ggplot2::unit(x = 0.075, units = 'inches')
        
      ), 
      
      show.legend = FALSE,
      linewidth = 1/2
      
    )
    
  }
  
  if(nrow(annotation_df_list$plus) > 0 | 
     nrow(annotation_df_list$minus) > 0){
    
    p <- p + 
      ggplot2::scale_color_manual(values = c('darkblue', 'darkgreen'))
    
  }
  
  
  
  p <- p + 
    ggplot2::xlim(start, end) +
    #ggplot2::scale_x_continuous(breaks = points, labels = pointlabels) +
    ggplot2::xlab(xlabel) + 
    ggplot2::ylab(ylabel) +
    ggplot2::geom_vline(xintercept = points[2], linetype = 2, color = 'red', size = 1) +
    ggplot2::facet_grid(condition~., scales = 'free_y') +
    ggplot2::ggtitle(binplottitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0)) +
    ggplot2::theme(panel.spacing = grid::unit(0, 'lines')) + 
    
    ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face),
                   #plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                   axis.title = ggplot2::element_text(size = titlesize, face = face),
                   axis.text = ggplot2::element_text(size = textsize, face = face),
                   strip.text = ggplot2::element_text(size = titlesize, face = face), 
                   
                   plot.margin = ggplot2::unit(rep(0.5, 4), 'cm'))
  
  print(p)
  
  
  
  
  #expandreads_totals <- expandplotdat$data
  
  p <- ggplot2::ggplot(expandreads_totals, mapping = ggplot2::aes(x = xcoord, y = reads))
  
  p <- p + ggplot2::geom_point(size = 2, color = scales::hue_pal()(1))
  
  ylims <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$continuous_range
  
  xlabel <- paste0(genedat$gene_id, ' ', chr, 
                   ':', expandannotation_df_list$labels$start, '-', 
                   expandannotation_df_list$labels$end, 'bp')
  
  if(expandannotation_df_list$geneplotlabel == 'Truncated'){
    
    xlabel <- paste0(genedat$gene_id, ' (Truncated) ', chr, 
                     ':', expandannotation_df_list$labels$start, '-', 
                     expandannotation_df_list$labels$end, 'bp')
    
  }
  
  
  
  p <- p + 
    
    # exons
    ggplot2::geom_segment(
      
      data = expandannotation_df_list$exons,
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        
        xend = end,
        #xend = expandplotend, 
        
        yend = ylims[1],
        color = strand
        
      ), 
      
      show.legend = FALSE,
      linewidth = 6
      
    ) + 
    
    # gene body
    ggplot2::geom_segment(
      
      data = expandannotation_df_list$labels,
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        xend = end,
        yend = ylims[1],
        color = strand
        
      ), 
      
      show.legend = FALSE,
      linewidth = 1
      
    )
  
  if(nrow(expandannotation_df_list$plus) > 0){
    # forward strand arrows
    p <- p + ggplot2::geom_segment(
      
      data = expandannotation_df_list$plus,
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        xend = end,
        yend = ylims[1],
        color = strand
        
      ), 
      
      arrow = ggplot2::arrow(
        
        ends = 'first',
        type = 'open',
        angle = 135,
        length = ggplot2::unit(x = 0.075, units = 'inches')
        
      ), 
      
      show.legend = FALSE,
      linewidth = 1/2
      
    )
    
  }  
  
  if(nrow(expandannotation_df_list$minus) > 0){
    
    # reverse strand arrows
    p <- p + ggplot2::geom_segment(
      
      data = expandannotation_df_list$minus, 
      
      mapping = ggplot2::aes(
        
        x = start,
        y = ylims[1],
        xend = end,
        yend = ylims[1],
        color = strand
        
      ), 
      
      arrow = ggplot2::arrow(
        
        ends = 'first',
        type = 'open',
        angle = 45,
        length = ggplot2::unit(x = 0.075, units = 'inches')
        
      ), 
      
      show.legend = FALSE,
      linewidth = 1/2
      
    )
    
  }
  
  if(nrow(expandannotation_df_list$plus) > 0 | 
     nrow(expandannotation_df_list$minus) > 0){
    
    p <- p + 
      ggplot2::scale_color_manual(values = c('darkblue', 'darkgreen'))
    
  }
  
  p <- p + 
    ggplot2::xlim(expandplotstart, expandplotend) +
    #ggplot2::scale_x_continuous(breaks = expandpoints, labels = expandpointlabels) +
    ggplot2::xlab(xlabel) + 
    ggplot2::ylab(ylabel) +
    ggplot2::geom_vline(xintercept = expandpoints[2], linetype = 1, color = 'red', size = 1) +
    ggplot2::facet_grid(condition~., scales = 'free_y') +
    ggplot2::ggtitle(expandplottitle) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0)) +
    ggplot2::theme(panel.spacing = grid::unit(0, 'lines')) + 
    
    ggplot2::theme(plot.title = ggplot2::element_text(size = titlesize, face = face),
                   #plot.subtitle = ggplot2::element_text(size = textsize, face = face),
                   axis.title = ggplot2::element_text(size = titlesize, face = face),
                   axis.text = ggplot2::element_text(size = textsize, face = face),
                   strip.text = ggplot2::element_text(size = titlesize, face = face), 
                   
                   plot.margin = ggplot2::unit(rep(0.5, 4), 'cm'))
  
  print(p)
  
  
  
}


