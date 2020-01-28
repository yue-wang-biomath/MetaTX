# MetaTX
## Introduction
The MetaTX is aimed for visualizing distribution of RNA-related genomic features in the transcriptome. 
To install exomePeak2 from Github, use the following codes.
```{r cars}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("SummarizedExperiment","cqn","Rsamtools",
                       "GenomicAlignments","GenomicRanges","GenomicFeatures",
                       "DESeq2","ggplot2","mclust",
                       "genefilter","BSgenome","BiocParallel",
                       "IRanges","S4Vectors","quantreg",
                       "reshape2","rtracklayer","apeglm","RMariaDB"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("yue-wang-biomath/MetaTX")
```
![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure1.png)
![image](https://github.com/yue-wang-biomath/MetaTX/blob/master/figure2.png)
# References
Chen, K. et al. (2014) High-resolution N6-methyladenosine (m6A) map using photo-crosslinking-assisted m6A sequencing. Angew. Chem. Int. Ed., 54, 1587-1590. (https://doi.org/10.1002/anie.201410647)

Schwartz, S. et al. (2014) Perturbation of m6A writers reveals two distinct classes of mRNA methylation at internal and 5'sites. Cell Reports, 8(1), 284-296. (https://doi.org/10.1016/j.celrep.2014.05.048)
