# MetaTX
## Introduction
The MetaTX is aimed for plotting the transcriptomic distribution of RNA-related genomic features.


## Quick Start with MetaTX

To install MetaTX from Github, please use the following codes.

```{r introduction}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicAlignments","GenomicRanges","GenomicFeatures", "ggplot2",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yue-wang-biomath/MetaTX")
```


## Visualization of the transcriptomic distribution 

It requires basic information of target feature set, involving the genomic locations, seqnames and strand types of each feature. The input feature set is required to be provided as a GRanges object. 

The ```metaTXplot``` function enables the visualization of RNA-related genomic features. First, the TxDb object need to be downloaded.

```
txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene

```

Load example dataset provided in the MetaTX package.

```
data("m6A_methyl_1")

```

Please see the following example, which will read 1000 m6A methylation sites from the file m6A_methyl_1 into R and sketch the distribution of these features along mRNA.

```
metaTXplot(m6A_methyl_1[1:1000], 
           txdb,  
           num_bin = 10,
           includeNeighborDNA = FALSE, 
           comparison = FALSE,
           lambda = 2,
           adjust = 0.15,
           title  = 'Distribution on mRNA',
           legend = 'PA-seq')            
``` 

![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure1.png)


The```metaTXplot``` function also enables the comparison between estimates via the MetaTX and the direct estimation method.

```
metaTXplot(m6A_methyl_1[1:1000], 
           txdb,  
           num_bin = 10,
           includeNeighborDNA = FALSE, 
           comparison = TRUE,
           lambda = 2,
           adjust = 0.15,
           title  = 'Distribution on mRNA',
           legend = 'PA-seq')            
``` 

![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure2.png)


## Resolving ambiguity problem

The package also provides an ```isoformProb``` function that can return the probabilities of a particular feature being located on different isoforms. 

```
isoformProb(features = m6A_methyl_1[1:1000],
            txdb,             
            num_bin = 10,
            includeNeighborDNA = TRUE)
```

Then probabilities of a particular feature being located on different isoforms (the last column) can be returned.

```
      index_trans index_methyl seqnames methyl_po strand trans_ID isoform_prob
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

Schwartz, S. et al. (2014) Perturbation of m6A writers reveals two distinct classes of mRNA methylation at internal and 5'sites. Cell Reports, 8(1), 284-296. (https://doi.org/10.1016/j.celrep.2014.05.048)
