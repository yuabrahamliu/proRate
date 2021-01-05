# Tutorial for R package proRate

### Yu Liu
### 1/3/2021

## Introduction

The dynamics of transcriptional elongation can influence various post-transcriptional processes, such as splicing, polyadenylation, and nuclear export, etc *(Bentley, 2014; Lei, et al., 2001; Wallace and Beggs, 2017)*. To quantify the elongation rate, a typical method is to treat cells with drugs able to inhibit Polymerase II (Pol II) from entering the genes and starting transcription, such as DRB (5,6-dichloro-1-β-d-ribofuranosylbenzimidazole), and then track Pol II using Pro-seq or Gro-seq *(Hou, et al., 2019)*. As the DRB dependent blocking of Pol II entry persists, a read blank region will form on gene body because few Pol II enter into the DNA template. Downstream of this blank region (Pol II depleted state) is the intact reads region (Pol II occupied state) formed by the unaffected Pol II having entered the gene before the drug blocking. Hence, the length of the blank region is considered as the distance that Pol II should have gone through during the blocking period, and the corresponding Pol II transcription rate can be obtained based on this distance and the blocking time.

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

All the 4 bam files are from a previous work studying the influence of the elongation factor Paf1 on RNA polymerase II transcription rate *(Hou, et al., 2019)*. Shortly, gene editing was used to produce mouse myoblast C2C12 cells with the factor Paf1 **conditionally** knocked out (KO cells), while some wild C2C12 cells with intact Paf1 expression were used as controls (WT cells). Then, both the KO cells and WT cells were divided into 2 groups, one of them was treated with the drug DRB to block RNA polymerase II from entering the gene templates and starting transcription for 15 min (15 min group), while the other was not treated with DRB (0 min group). Hence, the 4 bam files represent 4 experimental conditions. They are WT cells + 0 min DRB (wt0), WT cells + 15 min DRB (wt15), KO cells + 0 min DRB (ko0), and KO cells + 15 min DRB (ko15). After the DRB treatment step, Pro-seq was used to sequence nascent RNA for the C2C12 cells in these 4 conditions and the bam files were transformed from these Pro-seq data. To simplify the analysis in this tutorial, only reads from mouse chromosome X are included in the 4 bam files. 

The txt file that will be used in this tutorial records the coordinate information of 47 genes in mouse chromosome X, including their chromosomes (column "chr" in this txt file), their gene start coordinates (column "start" in this txt file), gene end coordinates (column "end"), strand +/- information (column "strand"), gene symbols (column "gene_id"), and corresponding mRNA accessions (column "id"). This file will provide these 47 genes as the ones need to be analyzed to *proRate*. It should be noted that the gene coordinates in this txt file are 0-based coordinates, which follows the rule of *Python*, but the functions in *proRate* will automatically convert them to 1-based coordinates internally, which follows the rule of *R*.

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

Hence, the transcription rate for the genes in WT C2C12 cells can be calculated using the following command.

```
wtrates <- calrate(time1file = wt0file, time2file = wt15file, 
                   time = 15, strandmethod = 1, 
                   targetfile = NULL, genomename = 'mm10', 
                   genelencutoff = 40000, fpkmcutoff = 1, 
                   savegenenames = c('Tbl1x', 'Gk'), 
                   plotgenenames = TRUE)
```

The result `wtrates` is a list and its most important element is a data.frame named "report". It contains the distance that Pol II passes (i.e. the length of the Pol II depleted region) on each gene template, which is inferred using the least sum of square (LSS) method of *proRate*. Additionally, the Pol II travel time (15 min here for all genes), the final rate calculated from these 2 values, the statistical significance of the results, as well as many other information are also included in this data.frame.

For the 2 genes Tbl1x and Gk, whose symbols have been transferred to the parameter `savegenenames`, their normalized reads values for each bp are saved in another 2 list elements named "genereads1" and "genereads2", which correspond to the values in `time1file` and `time2file` respectively. These values are also plotted by `calrate`, so the gene reads distributions in the 2 time files can be seen directly, as shown below.

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/tbl1xtracks_catch.PNG)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/Gktracks_catch.PNG)

For each gene, to infer its transition point between Pol II depleted region and occupied region, *proRate* first divides the gene body into 40 bins and calculates the reads ratio between `time2file` and `time1file` for each bin, and identifies the transition point on this bin level using the LSS method. The plots for the bin ratios for genes Tbl1x and Gk are also generated by `calrate`, as shown below.

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/tbl1xbins.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/Gkbins.png)

The transition point at this stage is actually a bin including hundreds of bases. To get a more precise position, this bin and its downstream one will be expanded to a single base resolution. The LSS method will be used on them to find the base with the smallest sum of square. This point is the final transition point. The single base level results on the 2 expanded bins are also plotted by `calrate`, for Tbl1x and Gk respectively.

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/tbl1xextends.png)

![](https://github.com/yuabrahamliu/proRate/blob/master/vignettes/Gkextends.png)

## Intra-rate statistics

After the rate inference for WT C2C12 cells, the transcriptional rates of different genes can be compared within the WT cells using the function `intracompare`. For its parameter `inferres`, transfer the result list from `calrate` directly to it, while for the parameter `quantilenum`, set it as 2 here, which means the genes will be divided into 2 quantile groups based on their transcriptional rates.

intracompareres <- intracompare(inferres = wtrates, quantilenum = 2)

The result `intracompareres` is a list with 2 elements. One is a data.frame containing the absolute rates (bp/min) of the genes, while the other records the relative rates (portion of the whole gene body/min). The genes with a greater absolute or relative rates are attributed as quantile2 genes, while the smaller ones are quantile1 genes.
