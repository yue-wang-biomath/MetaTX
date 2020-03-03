# MetaTX
## Introduction
The MetaTX is aimed for plotting the transcriptomic distribution of RNA-related genomic features.


## 1. Quick Start with MetaTX

To install MetaTX from Github, please use the following codes.

```{r introduction}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicAlignments","GenomicRanges","GenomicFeatures", "ggplot2",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"))
               
library(c("GenomicAlignments","GenomicRanges","GenomicFeatures", "ggplot2",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yue-wang-biomath/MetaTX")
library('MetaTX')
```


## 2. Data preprocessing 

First, the TxDb object need to be downloaded.

```
txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene
```
It requires basic information of target feature set, involving the genomic locations, seqnames and strand types of each feature. The input feature set is required to be provided as a GRanges object.

In the MetaTX package, we provide three example feature sets stored in the file ```m6A_methyl_1.rda```, ```m6A_methyl_2.rda``` and ```m6A_methyl_3.rda```. They are m6A datasets derived from different high-throughput sequencing approaches, including an miCLIP-seq dataset (Linder, et al., 2015; Olarerin-George and Jaffrey, 2017), a PA-m6A-seq dataset (Chen, et al., 2015) and an m6A-seq dataset (Schwartz, et al., 2014). 

We will use the example ```m6A_methyl_2.rda``` to illustrate how to use MetaTX sketching feature distribution. Users can also use other RNA-related genomic feature datasets.

Load the example dataset provided in the MetaTX package. 
```
data("m6A_methyl_2")
```

Please see the following example, which will read m6A methylation sites from the file ```m6A_methyl_2.rda``` into R and map these features to an mRNA model. 


```
remap_results_m6A <- remapCoord(features = m6A_methyl_2, txdb = txdb, num_bin = 10, includeNeighborDNA = TRUE) 
``` 

We provide this result and store it in the file ```remap_results_m6A.rda```.

## 3. Visualization of the transcriptomic distribution 

Use previously generated results or provided file ```remap_results_m6A.rda```.

```
data("remap_results_m6A")
```
The ```metaTXplot``` function enables the visualization of RNA-related genomic features.
```
txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene
metaTXplot(remap_results_m6A, txdb,  includeNeighborDNA = TRUE) 
``` 

![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure1.png)


The```metaTXplot``` function also enables the comparison between estimates via the MetaTX and the direct estimation method.

```
metaTXplot(remap_results_m6A, txdb, includeNeighborDNA = TRUE, comparison = TRUE)         
``` 

![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure2.png)


## 4. Resolving ambiguity problem

The package also provides an ```isoformProb``` function that can return the probabilities of a particular feature being located on different isoforms. 

```
data("remap_results_m6A")
txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene
isoform_probs <- isoformProb(remap_results_m6A, includeNeighborDNA = TRUE)

```

The probabilities of a particular feature being located on different isoforms (the last column) can be returned.

```
      index_trans index_methyl seqnames methyl_pos strand trans_ID isoform_prob
1             1            1    chr19     581474      +    65776 2.022115e-01
2             2            1    chr19     581474      +    65777 2.022115e-01
3             3            1    chr19     581474      +    65778 2.030733e-01
4             4            1    chr19     581474      +    65779 1.894304e-01
5             5            1    chr19     581474      +    65780 2.030733e-01
6             6            2     chr6  122744806      +    25251 5.000000e-01
7             7            2     chr6  122744806      +    25252 5.000000e-01
8             8            3     chr3  195250580      -    17278 1.000000e+00
```


# References
Chen, K. et al. (2014) High-resolution N6-methyladenosine (m6A) map using photo-crosslinking-assisted m6A sequencing. Angew. Chem. Int. Ed., 54, 1587-1590. (https://doi.org/10.1002/anie.201410647)

Linder,B. et al. (2015) Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. Nat. Methods, 12, 767–772. (https://doi.org/10.1038/nmeth.3453)

Olarerin-George, A.O. and Jaffrey, S.R. MetaPlotR: a Perl/R pipeline for plotting metagenes of nucleotide modifications and other transcriptomic sites. Bioinformatics 2017;33(10):1563-1564 (https://doi.org/10.1093/bioinformatics/btx002)

Schwartz, S. et al. (2014) Perturbation of m6A writers reveals two distinct classes of mRNA methylation at internal and 5'sites. Cell Reports, 8(1), 284-296. (https://doi.org/10.1016/j.celrep.2014.05.048)
