\name{fiveUTR_SIM}
\alias{fiveUTR_SIM}
\docType{data}
\title{
Sample data for checking CDS region
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("fiveUTR_SIM")}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 7 slots
  ..@ seqnames       :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 32 levels "chr1","chr1_gl000191_random",..: 6 18 12 3 17 30 10 13 8 17 ...
  .. .. ..@ lengths        : int [1:8349] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:8800] 49270717 57524031 50528607 99258051 47517368 134955272 21825373 103348864 77197699 32754761 ...
  .. .. ..@ width          : int [1:8800] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1 2 1 2 1 2 1 2 1 2 ...
  .. .. ..@ lengths        : int [1:4224] 1 5 1 2 1 1 5 3 1 1 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomeInfoDb"] with 4 slots
  .. .. ..@ seqnames   : chr [1:32] "chr1" "chr1_gl000191_random" "chr10" "chr11" ...
  .. .. ..@ seqlengths : int [1:32] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ is_circular: logi [1:32] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:32] NA NA NA NA ...
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
data(fiveUTR_SIM)
## maybe str(fiveUTR_SIM) ; plot(fiveUTR_SIM) ...
}
\keyword{datasets}
