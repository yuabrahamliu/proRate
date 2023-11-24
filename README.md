# Tutorial for R package proRate

### 
### 1/22/2021

## Introduction

The dynamics of transcriptional elongation can influence various post-transcriptional processes, such as splicing, polyadenylation, and nuclear export, etc *(Bentley, 2014; Lei, et al., 2001; Wallace and Beggs, 2017)*. To quantify the elongation rate, a typical method is to treat cells with drugs able to inhibit Polymerase II (Pol II) from entering the gene bodies and starting transcription, such as DRB (5,6-dichloro-1-β-d-ribofuranosylbenzimidazole), and then track Pol II using Pro-seq or Gro-seq *(Hou, et al., 2019)*. As the DRB dependent blocking of Pol II entry persists, a read blank region will form on gene body because few Pol II enter into it. Downstream of this blank region (Pol II depleted state) is the intact reads region (Pol II occupied state) formed by the unaffected Pol II having entered the gene body before the drug blocking. Hence, the length of the blank region is considered as the distance that Pol II should have gone through during the blocking period, and the corresponding Pol II transcription rate can be obtained based on this distance and the blocking time.

To get the distance, the transition point between the Pol II depleted state and occupied state must be identified. Traditionally, a 2-state hidden Markov model (HMM) is used to infer this elongation rate from Pro-seq or Gro-seq data *(Hou, et al., 2019)*. Here, we developed the R package *proRate*, which uses a novel least sum of square (LSS) method to calculate the rate and gets rid of the data distribution limitation of HMM model, as well as its complex parameter estimation steps. In addition, this package also offers other functions frequently used when study transcription dynamics, such as metagene plotting, pausing index calculation, differential gene identification, etc.

## Package installation

Code and documentation of *proRate* are freely available at <https://github.com/yuabrahamliu/proRate>, or <https://github.com/FADHLyemen/proRate>.

The following commands can be used to install this R package. 

```
library(devtools)

install_github('yuabrahamliu/proRate')
```

## Data preparation

To demonstrate the functions of *proRate*, this tutorial uses some data that accompany with the package, including 4 bam files and a txt file. 

All the 4 bam files are from a previous work studying the influence of the elongation factor Paf1 on RNA polymerase II transcription rate *(Hou, et al., 2019)*. Shortly, gene editing was used to produce mouse myoblast C2C12 cells with the factor Paf1 **conditionally** knocked out (KO cells), while some wild C2C12 cells with intact Paf1 expression were used as controls (WT cells). Then, both the KO cells and WT cells were divided into 2 groups, one of them was treated with the drug DRB to block RNA polymerase II from entering the gene bodies for 15 min (15 min group), while the other was not treated with DRB (0 min group). Hence, the 4 bam files represent 4 experimental conditions. They are WT cells + 0 min DRB (wt0), WT cells + 15 min DRB (wt15), KO cells + 0 min DRB (ko0), and KO cells + 15 min DRB (ko15). After the DRB treatment step, Pro-seq was used to sequence nascent RNA for the C2C12 cells in these 4 conditions and the bam files were transformed from these Pro-seq data. To simplify the analysis in this tutorial, only reads from mouse chromosome 19, which is the smallest autosome in mouse genome,  are included in the 4 bam files. 

The txt file that will be used in this tutorial records the coordinate information of 31 genes in mouse chromosome 19, including their chromosomes (column "chr" in this txt file), their gene start coordinates (column "start" in this txt file), gene end coordinates (column "end"), strand +/- information (column "strand"), gene symbols (column "gene_id"), and corresponding mRNA accessions (column "id"). This file will provide these 31 genes as the ones need to be analyzed to *proRate*. It should be noted that the gene coordinates in this txt file are 0-based coordinates, which follows the rule of *Python*, but the functions in *proRate* will automatically convert them to 1-based coordinates internally, which follows the rule of *R*.

Next, attach *proRate* to the R session and get the directories of its accompanied example files above.

```
library(proRate)

wt0file <- system.file('extdata', 'wt0.bam', package = 'proRate')
wt15file <- system.file('extdata', 'wt15.bam', package = 'proRate')
ko0file <- system.file('extdata', 'ko0.bam', package = 'proRate')
ko15file <- system.file('extdata', 'ko15.bam', package = 'proRate')

targetfile <- system.file('extdata', 'targetgenes.txt', package = 'proRate')
```

## Rate inference

We first use the function `calrate` to calculate the transcription elongation rate for genes in **WT** C2C12 cells. 

The most 2 important parameters are `time1file` and `time2file`. `time1file` corresponds to the Pro-seq bam file used as reference, i.e., with a DRB treatment time of 0 min. Here, the directory of wt0 file taken above should be transferred. While `time2file` corresponds to the bam file used as experimental group with a DRB treatment for more than 0 min. Here, the directory of wt15 file should be transferred. The time difference between these 2 files is 15 min, so for another parameter `time`, its value is 15. In addition, because the Pro-seq sequencing libraries of these bam files were prepared using *directional ligation method* *(Hou, et al., 2019)*, the parameter `strandmethod` should be set as 1, while in other cases, if *dUTP method* is used to construct libraries, it should be set as 2.

To define the target genes whose transcription rate need to be calculated, the `targetfile` directory taken above can be transferred to the parameter `targetfile` here. However, another way to define the genes is to set this parameter as NULL, but set the parameter `genomename` as "mm10", which means all the genes in the mouse "mm10" genome should be analyzed. In addition, the parameters `genelencutoff` and `fpkmcutoff`, which are set as 40000 and 1 respectively as an example, further require only the genes with a length greater than 40000bp, and an FPKM expression value greater than 1 in `time1file` will be analyzed.

Among the target genes defined, if there are some whose reads values need to be checked in detail, a vector containing the symbols of these genes can be transferred to the parameter `savegenenames`, and then their normalized reads values for each bp along the gene body will be returned with the final results. In addition, the `plotgenenames` parameter can be set as TRUE to futher plot and show the Pol II depleted region and occupied region of these genes.

Because bam file processing is time-consuming, the bam file relevant functions in *proRate*, such as `calrate` and `mcalrate`, provide a parallelization option and their parameter `threads` can be set to a multiple thread number to save the computing time. Its default value is 1.

Hence, the transcription rate for the genes in WT C2C12 cells can be calculated using the following command.

```
wtrates <- calrate(time1file = wt0file, time2file = wt15file, 
                   time = 15, strandmethod = 1, 
                   targetfile = NULL, genomename = 'mm10', 
                   genelencutoff = 40000, fpkmcutoff = 1, 
                   savegenenames = c('Mamdc2', 'Cpt1a'), 
                   plotgenenames = TRUE, 
				   threads = 1)
```

The result `wtrates` is a list and its most important element is a data.frame named "report". It contains the distance that Pol II passes (i.e. the length of the Pol II depleted region) on each gene, which is inferred using the least sum of square (LSS) method of *proRate*. Additionally, the Pol II travel time (15 min here for all genes), the final rate calculated from these 2 values, the statistical significance of the results, as well as many other information are also included in this data.frame.

For the 2 genes Mamdc2 and Cpt1a, whose symbols have been transferred to the parameter `savegenenames`, their normalized reads values for each bp are saved in another 2 list elements named "genereads1" and "genereads2", which correspond to the values in `time1file` and `time2file` respectively. These values are also plotted by `calrate`, so the gene reads distributions in the 2 time files can be seen directly, as shown below.

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/mamdc2tracks_catch.PNG)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/cpt1atracks_catch.PNG)

For each gene, to infer its transition point between Pol II depleted region and occupied region, *proRate* first divides the gene body into 40 bins and calculates the reads ratio between `time2file` and `time1file` for each bin, and identifies the transition point on this bin level using the LSS method. The plots for the bin ratios for genes Mamdc2 and Cpt1a are also generated by `calrate`, as shown below.

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/mamdc2bins.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/cpt1abins.png)

The transition point at this stage is actually a bin including hundreds of bases. To get a more precise position, this bin and its downstream one will be expanded to a single base resolution. The LSS method will be used on them to find the base with the smallest sum of square. This point is the final transition point. The single base level results on the 2 expanded bins are also plotted by `calrate`, for Mamdc2 and Cpt1a respectively.

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/mamdc2extends.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/cpt1aextends.png)

## Intra-rate statistics

After the rate inference for WT C2C12 cells, the transcriptional rates of different genes can be compared within the WT cells using the function `intracompare`. For its parameter `inferres`, transfer the result list from `calrate` directly to it, while for the parameter `quantilenum`, set it as 2 here, which means the genes will be divided into 2 quantile groups based on their transcriptional rates.

```
intracompareres <- intracompare(inferres = wtrates, quantilenum = 2)
```

The result `intracompareres` is a list with 2 elements. One is a data.frame containing the absolute rates (bp/min) of the genes, while the other records the relative rates (portion of the whole gene body/min). The genes with a greater absolute or relative rates are attributed as quantile2 genes, while the smaller ones are quantile1 genes.

## Gene structure

Because gene structure can influence transcriptional dynamics, *proRate* provides some functions working on this aspect. For example, `getgc` and `getexon` can calculate the GC contents and exon coverages for specific gene regions.

```
intraabsresgc <- getgc(genomename = 'mm10', 
                       genenames = intracompareres$absolute_rate$gene_id)

intraabsresexon <- getexon(genomename = 'mm10', 
                           genenames = intracompareres$absolute_rate$gene_id)
```

The parameter `genenames` defines which genes should be analyzed. 

To compare the GC content and exon coverage difference between the genes with an absolute rate on quantile2 level and that on quantile1 level, a function is written manually here and based on the results from `getgc` and `getexon`, it shows that the genes with a higher absolute rate (quantile2 genes) have a lower level in both GC content and exon coverage.

```
comparequantiles <- function(structureres, 
                             quantileres){ 
  
  library(ggplot2)
  
  dat <- merge(x = quantileres, y = structureres, by = c('gene_id'))
  dat <- dat[match(quantileres$gene_id, dat$gene_id),]
  dat <- dat[,c(1, 2, 3, ncol(dat))]
  
  yname <- names(dat)[ncol(dat)]
  
  names(dat) <- c('genesym', 'rate', 'quantile', 'metric')
  
  pval <- wilcox.test(subset(dat, quantile == 'quantile2')$metric, 
                      subset(dat, quantile == 'quantile1')$metric)$p.val
  
  p <- ggplot(data = dat, mapping = aes(x = quantile, y = metric, 
                                        fill = quantile))
  p <- p + geom_boxplot() + 
    ylab(yname) + 
    ggtitle(paste0(yname, ' (Wilcox p-val = ', signif(pval, 3), 
                   ')')) + 
    scale_fill_discrete(guide = FALSE) + 
    theme_bw()
  
  print(p)
  
  }

comparequantiles(structureres = intraabsresgc, 
                 quantileres = intracompareres$absolute_rate)
```

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/intragccomp_chr19.png)

```
comparequantiles(structureres = intraabsresexon, 
                 quantileres = intracompareres$absolute_rate)
```

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/intraexoncomp_chr19.png)

To check whether there is any k-mer difference between the gene body sequences of quantile2 genes and quantile1 genes, the function `getkmer` can be used. Its parameter `k` defines the length of the k-mer sequences need to be analyzed, and the parameters `genes1` and `genes2` require the gene symbols of the 2 sets need to be compared. Here, `k` is set as 6, and `genes1` and `genes2` are the quantile2 and quantile1 gene symbols respectively. The command is shown below.

```
intraabsreskmer <- getkmer(genomename = 'mm10', 
                           k = 6, 
                           genes1 = subset(intracompareres$absolute_rate, 
                                           absq == 'quantile2')$gene_id, 
                           genes2 = subset(intracompareres$absolute_rate, 
                                           absq == 'quantile1')$gene_id, 
                           feature = 'genebody')
```

The result `intraabsreskmer` is a data.frame containing all the 6-mers analyzed and their corresponding values on their frequency ratio between quantile2 and quantile1 gene sets, enrichment p-values, and adjusted p-values. 

A plotting function is written here to display the 6-mer ratio difference between the 2 gene sets.

```
plotkmer <- function(structureres = intraabsreskmer){
  
  library(ggplot2)
  
  dat <- structureres
  dat$order <- as.numeric(row.names(dat))
  dat$label <- dat$kmer
  dat$label[dat$order > 5] <- ''
  dat$color <- 'blue'
  dat$color[dat$order <= 5] <- 'red'
  
  k <- log(nrow(dat), base = 4)
  title <- paste0(nrow(dat), ' sorted ', k, 'mers')
  
  p <- ggplot(data = dat, mapping = aes(x = order, y = ratio, 
                                        label = label))
  p <- p + geom_point(color = dat$color) + 

    ggtitle(title) + 
    geom_text(hjust = 0, vjust = 0, color = 'red') + 
    theme_bw()
  
  print(p)

}

plotkmer(structureres = intraabsreskmer)
```

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/kmerplot_chr19.png)

## Multiple rates inference

The above function `calrate` can only infer gene elongation rates from one pair of bam files (wt15 vs wt0). If want to do the inference for multiple pairs of bam files at one time, the function `mcalrate` can be used.

For example, from 2 pairs of data (wt15 vs wt0, and ko15 vs ko0), the gene elongation rates in WT C2C12 cells and KO C2C12 cells can be obtained respectively.

The parameter `time1files` requires a vector with the reference bam file directories (0 min DRB treatment) as its elements. Here, it can be set as `c(wt0file, ko0file)`, while the parameter `time2files` requires the other vector containing the bam file directories for experimental groups. Here, it is `c(wt15file, ko15file)`. The paired bam files in `time1files` and `time2files` should have the same index in their vectors, so that the function can identify them as a pair. Then, another parameter `times` records the time difference for the WT group and KO group files. Because both of them are 15 min, its value is `c(15, 15)`

Although `wt0file` and `ko0file` are from different types of C2C12 cells, both of them are not treated with DRB, so no Pol II depleted region should appear on their gene bodies. If want to merge these 2 DRB 0 min files together and serve as a uniform reference for both the WT and KO group, the parameter `mergerefs` can be set as TRUE. Otherwise, set it as FALSE.

```
mrates <- mcalrate(time1files = c(wt0file, ko0file), 
                   time2files = c(wt15file, ko15file), 
                   time = c(15, 15), strandmethod = 1, 
                   mergerefs = FALSE, 
                   targetfile = NULL, genomename = 'mm10', 
                   genelencutoff = 40000, fpkmcutoff = 1, 
                   savegenenames = c('Mamdc2', 'Cpt1a'), 
                   plotgenenames = FALSE)
```

The result `mrates` is a list with 2 sub-lists. The first one is the result for the WT bam file pair, while the second one is for the KO bam file pair. For each of the sub-list, its structure is the same as that of the list generated by `calrate`, with a data.frame named "report" containing the main inference results, as well as 2 other elements recording the reads values for specific genes if some gene symbols have transferred to the parameter `savegenenames`.

## Inter-rate statistics

If want to check the elongation rate difference between the WT cells and the KO cells, the functions `orgintercompareinput` and `intercompare` can be used.

First, use `orgintercompareinput` to convert the result of `mcalrate` from the original list to a data.frame, and then, this data.frame can be transferred to `intercompare` via its parameter `inferresmat`, and the gene rate difference between the WT and KO cells, as suggested by the parameter `groupname`, can be calculated.

```
orgmrates <- orgintercompareinput(reslist = mrates)

intercompareres <- intercompare(inferresmat = orgmrates, 
                                groupnames = c('WT', 'KO'))
```

The result `intercompareres` is a data.frame showing the elongation rate in WT cells and KO cells for each gene detected, as well as the p-value, adjusted p-value, and log2FC(KO/WT) for their KO/WT difference. Because both the WT and KO group here only has one experimental replicate, this time the p-value is calculated using Fisher's test without considering the variance of the observation.

## Metagene plotting

For transcriptional dynamic study, metagene plotting is always needed to check the gene reads distribution on a metagene level. Hence, *proRate* contains the functions `metaplot` and `mmetaplot` to do this. If want to include the metagene curves of multiple experimental conditions, such as the WT0 and KO0 here, into one plot, `mmetaplot` should be used.

For its parameter `metafiles`, the bam file directories of different experimental conditions should be transferred as a vector, so here set it as `c(wt0file, ko0file)`. Correspondingly, the parameter `labels` should be set as `c('WT', 'KO')`, so that the final metagene curves can be labeled with these conditions in the plot.

The final metagene plot can focus on different gene regions. If want to focus on promoter, the parameter `tssradius` can be set using a single numeric value or a numeric vector. For example, if set it with the vector `c(1000, 500)`, then 2 metagene plots will be generated. One covers the promoter region from 1000bp upstream to 1000bp downstream of the TSS point, the other covers the promoter from 500bp upstream to 500bp downstream of TSS. The parameter `ttsradius` can also be set in this way, if want to focus on the gene tail region.

If need to generate a metagene plot on the whole gene body, the parameter `genebodylen` should be given a numeric value, for example 2000. Different from the promoter region and gene tail region, which can be set a uniform length for all genes and then summarized together to the metagene level, the length of the whole gene body varies largely for different genes. Hence, this will be adjusted according to the value of `genebodylen`. If a gene has a length greater than the defined value of 2000bp here, its gene body will be compressed, while if a gene is shorter than 2000bp, it will be extended, so finally all the genes will be scaled to 2000bp, and the metagene plot can be generated. The command is shown below.

```
metareslist <- mmetaplot(metafiles = c(wt0file, ko0file), 
                         labels = c('WT', 'KO'), strandmethod = 1, 
                         targetgenefile = NULL, 
                         genomename = 'mm10', 
                         genelencutoff = 40000, 
                         fpkmcutoff = 1, 
                         tssradius = c(1000, 500), ttsradius = c(1000), 
                         genebodylen = 2000,
                         savegenenames = c('Mamdc2', 'Cpt1a'), 
						 plotgenenames = TRUE)
```

Several metagene plots can be returned, focusing on the various gene regions defined, as shown below. Because 2 specific genes, Mamdc2 and Cpt1a, have been transferred to the parameter `savegenenames`, their reads distributions are also plotted and returned.

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/metatss1000_chr19.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/metatss500_chr19.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/metatts1000_chr19.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/metagenebody_chr19.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/metagenebody1000_chr19.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/metamamdc2_catch.PNG)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/metacpt1a_catch.PNG)

The list result `metareslist` records the FPKM values for various gene regions. For example, its slot `TSSFPKMstatses` contains a data.frame with each row representing a gene and each column representing a promoter region, on either the transcriptional forward (sense) or reverse (antisense) strand. The corresponding FPKM values are recorded in this data.frame.

The exact FPM values shown in the metagene plots are also included in `metareslist`. For examples, `metareslist$WT$TSSfwdmeans` is a vector with a length of 2001, covering the FPM values from 1000bp upstream of TSS, then to the single TSS point, and finally to 1000bp downstream of TSS, on the forward (sense) strand.

## Pause index

Pause index (or pausing index) is the ratio of transcription polymerase II signal density near a gene promoter to signal density within the gene body *(Williams, et al., 2015)*. A gene with a large pause index means most polymerase II molecules accumulate in its promoter region, while few of them enter into the downstream gene body. Hence, most polymerases pause in the promoter and cannot be released to complete the downstream transcription.

The functions `calpauseidx` and `mcalpauseidx` in *proRate* can be used to calculate it. The command of `mcalpauseidx` is shown below to get gene pause indeces for multiple conditions (WT0 and KO0) at one time. The parameter `bamfiles` receives a vector with the bam file directories for different conditions. Here, it is `c(wt0file, ko0file)`. The parameter `labels` indicates the names of these conditions as `c('WT', 'KO')`.

To define the length of promoter, `tssradius` is set as 1000, so that the region from 1000bp upstream to 1000bp downstream of TSS is defined as promoter, and the further downstream region until TTS is defined as within gene body. The polymerase II densities of these 2 parts will be used to calculate pause index for each gene.

```
pauselist <- mcalpauseidx(bamfiles = c(wt0file, ko0file), 
                          labels = c('WT', 'KO'), 
                          strandmethod = 1, 
                          genefile = NULL, 
                          genomename = 'mm10', 
                          genelencutoff = 40000, fpkmcutoff = 1, 
                          tssradius = 1000, 
                          threads = 1)
```

The result `pauselist` contains 2 sub-lists. One is for WT condition, the other is for KO condition, both of which further include a data.frame named "report" recording the gene pause indeces for the WT or KO condition. Each row is a gene, while the columns named "pausedens" and "elongdens" are the polymerase II densities for the promoter region and the downstream gene body region. The next column "pauseidx" is the final pause indeces calculated from them. Then, the columns "pval" and "padj" reflect the statistical significance of this result judged using Fisher's test. 

## References

Bentley, D.L. Coupling mRNA processing with transcription in time and space. Nature Reviews Genetics 2014;15(3):163-175.

Hou, L., et al. Paf1C regulates RNA polymerase II progression by modulating elongation rate. Proceedings of the National Academy of Sciences 2019;116(29):14583-14592.

Lei, E.P., Krebber, H. and Silver, P.A. Messenger RNAs are recruited for nuclear export during transcription. Genes & Development 2001;15(14):1771-1782.

Wallace, E.W.J. and Beggs, J.D. Extremely fast and incredibly close: cotranscriptional splicing in budding yeast. RNA 2017;23(5):601-610.

Williams, Lucy H., et al. Pausing of RNA Polymerase II Regulates Mammalian Developmental Potential through Control of Signaling Networks. Molecular Cell 2015;58(2):311-322.


