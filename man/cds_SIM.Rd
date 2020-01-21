\name{cds_SIM}
\alias{cds_SIM}
\docType{data}
\title{
AA}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("cds_SIM")}
\format{
  The format is:
Formal class 'GRanges' [package "GenomicRanges"] with 7 slots
  ..@ seqnames       :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 34 levels "chr1","chr1_gl000191_random",..: 30 15 5 26 24 5 33 25 9 34 ...
  .. .. ..@ lengths        : int [1:9891] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ ranges         :Formal class 'IRanges' [package "IRanges"] with 6 slots
  .. .. ..@ start          : int [1:10500] 158704321 43117943 15742399 4599620 1152750 70188277 117566836 2961268 71799478 25380631 ...
  .. .. ..@ width          : int [1:10500] 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ NAMES          : NULL
  .. .. ..@ elementType    : chr "ANY"
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ strand         :Formal class 'Rle' [package "S4Vectors"] with 4 slots
  .. .. ..@ values         : Factor w/ 3 levels "+","-","*": 1 2 1 2 1 2 1 2 1 2 ...
  .. .. ..@ lengths        : int [1:5208] 3 1 4 1 1 1 1 2 4 7 ...
  .. .. ..@ elementMetadata: NULL
  .. .. ..@ metadata       : list()
  ..@ seqinfo        :Formal class 'Seqinfo' [package "GenomeInfoDb"] with 4 slots
  .. .. ..@ seqnames   : chr [1:34] "chr1" "chr1_gl000191_random" "chr10" "chr11" ...
  .. .. ..@ seqlengths : int [1:34] NA NA NA NA NA NA NA NA NA NA ...
  .. .. ..@ is_circular: logi [1:34] NA NA NA NA NA NA ...
  .. .. ..@ genome     : chr [1:34] NA NA NA NA ...
  ..@ elementMetadata:Formal class 'DataFrame' [package "S4Vectors"] with 6 slots
  .. .. ..@ rownames       : NULL
  .. .. ..@ nrows          : int 10500
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
data(cds_SIM)
## maybe str(cds_SIM) ; plot(cds_SIM) ...
}
\keyword{datasets}
