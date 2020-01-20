\name{get_promoter_function}
\alias{get_promoter_function}
Function: Map features to the promoter region
\title{
Map features to the promoter region
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_promoter_function(width_promoter, num_bin_promoter, fiveUTR_align, fiveUTR_tx0, methyl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{width_promoter}{
%%     ~~Describe \code{width_promoter} here~~
}
  \item{num_bin_promoter}{
%%     ~~Describe \code{num_bin_promoter} here~~
}
  \item{fiveUTR_align}{
%%     ~~Describe \code{fiveUTR_align} here~~
}
  \item{fiveUTR_tx0}{
%%     ~~Describe \code{fiveUTR_tx0} here~~
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
function (width_promoter, num_bin_promoter, fiveUTR_align, fiveUTR_tx0, 
    methyl) 
{
    fiveUTR_align_non_exist <- fiveUTR_align[which(fiveUTR_align[, 
        "width"] == 0), ]
    fiveUTR_align_exist <- fiveUTR_align[which(fiveUTR_align[, 
        "width"] != 0), ]
    fiveUTR_align_pstv <- fiveUTR_align_exist[fiveUTR_align_exist[, 
        "strand"] == "+", ]
    fiveUTR_align_ngtv <- fiveUTR_align_exist[fiveUTR_align_exist[, 
        "strand"] == "-", ]
    fiveUTR_pstv <- fiveUTR_tx0[as.character(fiveUTR_align_pstv[, 
        "trans_ID"])]
    fiveUTR_ngtv <- fiveUTR_tx0[as.character(fiveUTR_align_ngtv[, 
        "trans_ID"])]
    promoter_pstv_start <- GRanges()
    promoter_ngtv_start <- GRanges()
    if (nrow(fiveUTR_align_pstv) != 0) {
        promoter_pstv_start <- GRanges(seqnames = fiveUTR_align_pstv[, 
            "seqnames"], ranges = IRanges(min(start(ranges(fiveUTR_pstv))), 
            width = 1), strand = "+")
    }
    if (nrow(fiveUTR_align_ngtv) != 0) {
        promoter_ngtv_start <- GRanges(seqnames = fiveUTR_align_ngtv[, 
            "seqnames"], ranges = IRanges(max(end(ranges(fiveUTR_ngtv))), 
            width = 1), strand = "-")
    }
    promoter_pstv_rgs <- flank(promoter_pstv_start, width_tail, 
        start = TRUE)
    promoter_ngtv_rgs <- flank(promoter_ngtv_start, width_tail, 
        start = TRUE)
    methyl_promoter_map_pstv <- data.frame(findOverlaps(methyl, 
        promoter_pstv_rgs))
    methyl_promoter_map_ngtv <- data.frame(findOverlaps(methyl, 
        promoter_ngtv_rgs))
    dist_from_start <- start(methyl[methyl_promoter_map_pstv[, 
        "queryHits"]]) - start(promoter_pstv_rgs)[methyl_promoter_map_pstv[, 
        "subjectHits"]]
    align_mtr_index <- ceiling(dist_from_start/width_promoter * 
        num_bin_promoter)
    align_mtr_pstv <- matrix(0, nrow(fiveUTR_align_pstv), num_bin_promoter)
    for (i in 1:nrow(align_mtr_pstv)) {
        align_mtr_pstv[methyl_promoter_map_pstv[, "subjectHits"][i], 
            align_mtr_index[i]] <- 1
    }
    promoter_align_pstv <- data.frame(fiveUTR_align_pstv[, c("index_trans", 
        "index_methyl", "seqnames", "methyl_pos", "strand", "trans_ID")], 
        width = replicate(nrow(align_mtr_pstv), width_promoter), 
        coordinate = align_mtr_pstv)
    dist_from_start <- end(promoter_ngtv_rgs)[methyl_promoter_map_ngtv[, 
        "subjectHits"]] - start(methyl[methyl_promoter_map_ngtv[, 
        "queryHits"]])
    align_mtr_index <- ceiling(dist_from_start/width_promoter * 
        num_bin_promoter)
    align_mtr_ngtv <- matrix(0, nrow(fiveUTR_align_ngtv), num_bin_promoter)
    for (i in 1:nrow(align_mtr_ngtv)) {
        align_mtr_ngtv[methyl_promoter_map_ngtv[, "subjectHits"][i], 
            align_mtr_index[i]] <- 1
    }
    promoter_align_ngtv <- data.frame(fiveUTR_align_ngtv[, c("index_trans", 
        "index_methyl", "seqnames", "methyl_pos", "strand", "trans_ID")], 
        width = replicate(nrow(align_mtr_ngtv), width_promoter), 
        coordinate = align_mtr_ngtv)
    promoter_align_non_exist <- data.frame(fiveUTR_align_non_exist[, 
        c("index_trans", "index_methyl", "seqnames", "methyl_pos", 
            "strand", "trans_ID", "width")], coordinate = matrix(0, 
        nrow(fiveUTR_align_non_exist), num_bin_promoter))
    promoter_align <- rbind(promoter_align_pstv, promoter_align_ngtv, 
        promoter_align_non_exist)
    promoter_align <- promoter_align[order(promoter_align[, "index_trans"]), 
        ]
    return(promoter_align)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
