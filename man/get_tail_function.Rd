\name{get_tail_function}
\alias{get_tail_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function: Map features to the tail region
}
\description{
 Map features to the tail region
}
\usage{
get_tail_function(width_tail, num_bin_tail, threeUTR_align, threeUTR_tx0, methyl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{width_tail}{
%%     ~~Describe \code{width_tail} here~~
}
  \item{num_bin_tail}{
%%     ~~Describe \code{num_bin_tail} here~~
}
  \item{threeUTR_align}{
%%     ~~Describe \code{threeUTR_align} here~~
}
  \item{threeUTR_tx0}{
%%     ~~Describe \code{threeUTR_tx0} here~~
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
function (width_tail, num_bin_tail, threeUTR_align, threeUTR_tx0, 
    methyl) 
{
    threeUTR_align_non_exist <- threeUTR_align[which(threeUTR_align[, 
        "width"] == 0), ]
    threeUTR_align_exist <- threeUTR_align[which(threeUTR_align[, 
        "width"] != 0), ]
    threeUTR_align_pstv <- threeUTR_align_exist[threeUTR_align_exist[, 
        "strand"] == "+", ]
    threeUTR_align_ngtv <- threeUTR_align_exist[threeUTR_align_exist[, 
        "strand"] == "-", ]
    threeUTR_pstv <- threeUTR_tx0[as.character(threeUTR_align_pstv[, 
        "trans_ID"])]
    threeUTR_ngtv <- threeUTR_tx0[as.character(threeUTR_align_ngtv[, 
        "trans_ID"])]
    tail_pstv_start <- GRanges()
    tail_ngtv_start <- GRanges()
    if (nrow(threeUTR_align_pstv) != 0) {
        tail_pstv_start <- GRanges(seqnames = threeUTR_align_pstv[, 
            "seqnames"], ranges = IRanges(max(end(ranges(threeUTR_pstv))), 
            width = 1), strand = "+")
    }
    if (nrow(threeUTR_align_ngtv) != 0) {
        tail_ngtv_start <- GRanges(seqnames = threeUTR_align_ngtv[, 
            "seqnames"], ranges = IRanges(min(start(ranges(threeUTR_ngtv))), 
            width = 1), strand = "-")
    }
    tail_pstv_rgs <- flank(tail_pstv_start, width_tail, start = FALSE)
    tail_ngtv_rgs <- flank(tail_ngtv_start, width_tail, start = FALSE)
    methyl_tail_map_pstv <- data.frame(findOverlaps(methyl, tail_pstv_rgs))
    methyl_tail_map_ngtv <- data.frame(findOverlaps(methyl, tail_ngtv_rgs))
    dist_from_start <- start(methyl[methyl_tail_map_pstv[, "queryHits"]]) - 
        start(tail_pstv_rgs)[methyl_tail_map_pstv[, "subjectHits"]]
    align_mtr_index <- ceiling(dist_from_start/width_tail * num_bin_tail)
    align_mtr_pstv <- matrix(0, nrow(threeUTR_align_pstv), num_bin_tail)
    for (i in 1:nrow(align_mtr_pstv)) {
        align_mtr_pstv[methyl_tail_map_pstv[, "subjectHits"][i], 
            align_mtr_index[i]] <- 1
    }
    tail_align_pstv <- data.frame(threeUTR_align_pstv[, c("index_trans", 
        "index_methyl", "seqnames", "methyl_pos", "strand", "trans_ID")], 
        width = replicate(nrow(align_mtr_pstv), width_tail), 
        coordinate = align_mtr_pstv)
    dist_from_start <- end(tail_ngtv_rgs)[methyl_tail_map_ngtv[, 
        "subjectHits"]] - start(methyl[methyl_tail_map_ngtv[, 
        "queryHits"]])
    align_mtr_index <- ceiling(dist_from_start/width_tail * num_bin_tail)
    align_mtr_ngtv <- matrix(0, nrow(threeUTR_align_ngtv), num_bin_tail)
    for (i in 1:nrow(align_mtr_ngtv)) {
        align_mtr_ngtv[methyl_tail_map_ngtv[, "subjectHits"][i], 
            align_mtr_index[i]] <- 1
    }
    tail_align_ngtv <- data.frame(threeUTR_align_ngtv[, c("index_trans", 
        "index_methyl", "seqnames", "methyl_pos", "strand", "trans_ID")], 
        width = replicate(nrow(align_mtr_ngtv), width_tail), 
        coordinate = align_mtr_ngtv)
    tail_align_non_exist <- data.frame(threeUTR_align_non_exist[, 
        c("index_trans", "index_methyl", "seqnames", "methyl_pos", 
            "strand", "trans_ID", "width")], coordinate = matrix(0, 
        nrow(threeUTR_align_non_exist), num_bin_tail))
    tail_align <- rbind(tail_align_pstv, tail_align_ngtv, tail_align_non_exist)
    tail_align <- tail_align[order(tail_align[, "index_trans"]), 
        ]
    return(tail_align)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
