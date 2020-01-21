\name{isoformProb}
\alias{isoformProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
isoformProb(features, txdb, num_bin = 10, includeNeighborDNA = FALSE)
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
    width_mtr <- remap_results[[2]]
    trans_info <- remap_results[[3]]
    num_bin_sum <- ncol(align_mtr)
    if (includeNeighborDNA) {
        weight_mtr <- replicate(num_bin_sum, rowSums(width_mtr))
        alpha <- get_correct_prob_function(num_bin_sum, align_mtr, 
            weight_mtr^lambda, trans_info)[[1]]
        prob_mtr <- get_correct_prob_function(num_bin_sum, align_mtr, 
            weight_mtr^lambda, trans_info)[[2]]
    }
    else {
        weight_start <- num_bin_sum/5 + 1
        weight_end <- num_bin_sum * 4/5
        weight_mtr <- replicate(num_bin_sum, rowSums(width_mtr[, 
            weight_start:weight_end]))
        alpha <- get_correct_prob_function(num_bin_sum, align_mtr, 
            weight_mtr^lambda, trans_info)[[1]]
        prob_mtr <- get_correct_prob_function(num_bin_sum, align_mtr, 
            weight_mtr^lambda, trans_info)[[2]]
    }
    prob <- as.numeric(rowSums(prob_mtr))
    prob_isoforms <- cbind(trans_info, isoform_prob = prob)
    return(prob_isoforms)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
