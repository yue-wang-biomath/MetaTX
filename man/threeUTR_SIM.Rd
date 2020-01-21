\name{threeUTR_SIM}
\alias{threeUTR_SIM}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("threeUTR_SIM")}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 7 slots
  ..@ seqnames       :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 33 levels "chr1","chr10",..: 32 3 18 2 12 6 18 30 1 27 ...
  .. .. ..@ lengths        : int [1:8437] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:8800] 150156635 113281102 160968954 101370279 10942579 54942255 156424006 116303794 169356440 2857307 ...
  .. .. ..@ width          : int [1:8800] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1 2 1 2 1 2 1 2 1 2 ...
  .. .. ..@ lengths        : int [1:4290] 1 1 1 1 1 1 4 6 1 1 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomeInfoDb"] with 4 slots
  .. .. ..@ seqnames   : chr [1:33] "chr1" "chr10" "chr11" "chr12" ...
  .. .. ..@ seqlengths : int [1:33] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ is_circular: logi [1:33] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:33] NA NA NA NA ...
  ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 8800
  .. .. ..@ listData       : Named list()
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ elementType    : chr "ANY"
  ..@ metadata       : list()
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(threeUTR_SIM)
## maybe str(threeUTR_SIM) ; plot(threeUTR_SIM) ...
}
\keyword{datasets}
