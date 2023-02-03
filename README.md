1.	Description  
AdaptRM has a multi-tasking framework for the synergetic learning of epitranscriptomes. It was enabled by the adaptive pooling and several standard convolution blocks. Trained for an integrated multi-task formulated by three case studies, the model can operate on both low-resolution and high-resolution datasets without further preprocessing the input primary sequence, and conduct tissue-specific and type-specific prediction according to users’ needs. 

2.	Requirements  
Before prediction, please make sure the following packages are installed in the Python environment:  
python = 3.10.0  
pytorch = 1.12.0  
numpy = 1.22.3  
pandas = 1.4.4  
argparse = 1.4.0  
scikit-learn = 1.1.1  
autograd = 1.5.0  
listmodel = 0.2.1  

3.	Input and output  
* The input object should be a list of [1, -1, 4] representing a query sequence.  
* The output would be a vector of values. Each element ranges from 0 to 1 and represents the probability of its corresponding task.  

  References  
Chen, K. et al. (2014) High-resolution N6-methyladenosine (m6A) map using photo-crosslinking-assisted m6A sequencing. Angew. Chem. Int. Ed., 54, 1587-1590. (https://doi.org/10.1002/anie.201410647)

Linder,B. et al. (2015) Single-nucleotide-resolution mapping of m6A and m6Am throughout the transcriptome. Nat. Methods, 12, 767–772. (https://doi.org/10.1038/nmeth.3453)

Olarerin-George, A.O. and Jaffrey, S.R. MetaPlotR: a Perl/R pipeline for plotting metagenes of nucleotide modifications and other transcriptomic sites. Bioinformatics 2017;33(10):1563-1564 (https://doi.org/10.1093/bioinformatics/btx002)

Schwartz, S. et al. (2014) Perturbation of m6A writers reveals two distinct classes of mRNA methylation at internal and 5'sites. Cell Reports, 8(1), 284-296. (https://doi.org/10.1016/j.celrep.2014.05.048)
