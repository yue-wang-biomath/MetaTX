\name{get_trans_function}
\alias{get_trans_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_trans_function(trans, methyl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trans}{
%%     ~~Describe \code{trans} here~~
}
  \item{methyl}{
%%     ~~Describe \code{methyl} here~~
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
function (trans, methyl) 
{
    sample_methyl_pos <- start(ranges(methyl))
    sample_methyl_seqnames <- data.frame(seqnames(methyl))[[1]]
    trans_metadata <- mcols(trans)
    index_trans <- 1:length(trans)
    index_methyl <- data.frame(trans)[, "xHits"]
    seqnames <- sample_methyl_seqnames[data.frame(trans)[, "xHits"]]
    methyl_pos <- sample_methyl_pos[data.frame(trans)[, "xHits"]]
    strand <- data.frame(trans)[, "strand"]
    trans_ID <- data.frame(trans)[, "transcriptsHits"]
    trans_info <- data.frame(index_trans = index_trans, index_methyl = index_methyl, 
        seqnames = seqnames, methyl_pos = methyl_pos, strand = strand, 
        trans_ID = trans_ID)
    return(trans_info)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
