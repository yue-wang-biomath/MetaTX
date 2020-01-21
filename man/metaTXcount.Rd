\name{metaTXcount}
\alias{metaTXcount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metaTXcount(features, txdb, num_bin = 10, includeNeighborDNA = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{features}{
%%     ~~Describe \code{features} here~~
}
  \item{txdb}{
%%     ~~Describe \code{txdb} here~~
}
  \item{num_bin}{
%%     ~~Describe \code{num_bin} here~~
}
  \item{includeNeighborDNA}{
%%     ~~Describe \code{includeNeighborDNA} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (features, txdb, num_bin = 10, includeNeighborDNA = FALSE) 
{
    remap_results <- remapCoord(methyl, txdb, num_bin, includeNeighborDNA)
    align_mtr <- remap_results[[1]]
    num_bin_sum <- ncol(align_mtr)
    feature_counts <- colSums(align_mtr)
    return(feature_counts)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
