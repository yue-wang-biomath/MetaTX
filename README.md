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
The```metaTXplot``` function enables the comparison between estimates via the MetaTX and the direct estimation method.
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

## Visualization of the transcriptomic distribution 

# References
Chen, K. et al. (2014) High-resolution N6-methyladenosine (m6A) map using photo-crosslinking-assisted m6A sequencing. Angew. Chem. Int. Ed., 54, 1587-1590. (https://doi.org/10.1002/anie.201410647)

Linder,B. et al. (2015) Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. Nat. Methods, 12, 767–772. (https://doi.org/10.1038/nmeth.3453)

Schwartz, S. et al. (2014) Perturbation of m6A writers reveals two distinct classes of mRNA methylation at internal and 5'sites. Cell Reports, 8(1), 284-296. (https://doi.org/10.1016/j.celrep.2014.05.048)
