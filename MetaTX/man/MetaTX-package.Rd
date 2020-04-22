\name{MetaTX-package}
\alias{MetaTX-package}
\alias{MetaTX}
\docType{package}
\title{
\packageTitle{MetaTX}
}
\description{
\packageDescription{MetaTX}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MetaTX}
\packageIndices{MetaTX}
}
\author{
\packageAuthor{MetaTX}

Maintainer: \packageMaintainer{MetaTX}
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
~~ documentation directory ~~
\keyword{ package }
\examples{
## load remap_results_m6A
data(remap_results_m6A)

## download the txdb object from UCSC
txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene

## plot
p1 <- metaTXplot(remap_results_m6A, txdb,  
           includeNeighborDNA = TRUE, 
           comparison = TRUE)
p1
}
