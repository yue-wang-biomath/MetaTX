\name{remapCoord}
\alias{remapCoord}
\title{
remapCoord
}

\description{

Remap genomic features represented by genome-based coordinate to an mRNA model. An alignment matrix, a width matrix and other associated annotations can be returned.

}

\usage{

remapCoord(features, txdb = NA, num_bin = 10, includeNeighborDNA = FALSE)

}

\arguments{

  \item{features}{
Genomic features, which should be a GRanges object.
}

  \item{txdb}{
A TxDb object, which should be a container for storing transcript annotations.
}

  \item{num_bin}{
The number of bins into which an mRNA component is divided.
}

  \item{includeNeighborDNA}{
Whether neighborhood DNA regions should be included.
}

}

\details{

This function returns a list contains three components, which are an alignment matrix, a width matrix, and a matrix of other associated annotations.

}

\value{

A list object. 

}

\examples{
txdb  <- TxDb.Hsapiens.UCSC.hg19.knownGene

## read genomic features:
data("m6A_methyl_1")

## remap genomic features to an mRNA model:
remap_results <- remapCoord(features = m6A_methyl_1, txdb = txdb, num_bin = 10)

## include neighborhood DNA regions:
remap_results_m6A <- remapCoord(features = m6A_methyl_1, txdb = txdb, num_bin = 10, includeNeighborDNA = TRUE)

}
