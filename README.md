# proRate

This program provides functions used to infer RNA Polymerase 2 (Pol2) elongation rate from Pro-seq/Gro-seq bam files using a least sum of square (LSS) method.

## Usage

To calculate Pol2 elongation rate, 2 files Pro-seq/Gro-seq bam files are needed. One is generated from sample treated with DRB (5,6-Dichloro-1-beta-D-ribofuranosylbenzimidazole, transcriptional elongation inhibitor) for a specific time period (e.g. 15min). The other is generated from reference sample without DRB treatment and deemed as a sample for the time point of 0min.

After uncompress the "elongHMM" package, set the parameters in the "config" file, and then put time1file, time2file and targetfile mentioned below in the "elongHMM/elongHMM_pipeline/data" directory. Then, in the elongHMM directory, simply type `make`, the program will run automatically. The result will be found in the "elongHMM/elongHMM_pipeline/result" directory.

## Paramters

* time1file - The reference sample Pro-seq/Gro-seq name, corresponding to 0min.
* time2file - The treatment sample Pro-seq/Gro-seq name, corresponding to a specific time point (e.g. 15min if DRB was used for 15min)
* time - The time difference between time1file and time2file, using min as its unit.
* samplepair - The name of the time1file/time2file sample pair. It will be used as the prefix of the final result file.
* targetfile - The genes whose transcriptional elongation rate need to be analyzed. The chrom, start, end, strand, sym and id information of the genes need to be specified in this file. Refer to the example targetfile named "uniq_target.txt".
* strandmethod - The strand specific method used when preparing the sequencing library, can be 1 for directional ligation method and 2 for dUTP method. If the sample is sequenced using a single strand method, just set it as "single".

