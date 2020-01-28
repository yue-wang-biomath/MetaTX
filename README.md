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

The ```metaTXplot``` function enables the visualization of RNA-related genomic features. Please see the following example, which will read 9000 m6A methylation sites from the file m6A_methyl_1 into R and sketch the distribution of these features along mRNA.


![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure1.png)
![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure2.png)
# References
Chen, K. et al. (2014) High-resolution N6-methyladenosine (m6A) map using photo-crosslinking-assisted m6A sequencing. Angew. Chem. Int. Ed., 54, 1587-1590. (https://doi.org/10.1002/anie.201410647)

Schwartz, S. et al. (2014) Perturbation of m6A writers reveals two distinct classes of mRNA methylation at internal and 5'sites. Cell Reports, 8(1), 284-296. (https://doi.org/10.1016/j.celrep.2014.05.048)
