\name{remapCoord}
\alias{remapCoord}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
remapCoord(features, txdb = NA, num_bin = 10, includeNeighborDNA = FALSE)
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
function (features, txdb = NA, num_bin = 10, includeNeighborDNA = FALSE) 
{
    cds_by_tx0 <- cdsBy(txdb, "tx")
    fiveUTR_tx0 <- fiveUTRsByTranscript(txdb, use.names = FALSE)
    threeUTR_tx0 <- threeUTRsByTranscript(txdb, use.names = FALSE)
    methyl <- GRanges(Rle(seqnames(features)), IRanges(as.numeric(end(features)), 
        width = 1))
    strand(methyl) <- Rle(strand(features))
    trans <- mapToTranscripts(methyl, txdb)
    trans_info <- get_trans_function(trans, methyl)
    num_bin_fiveUTR <- num_bin
    num_bin_cds <- num_bin
    num_bin_threeUTR <- num_bin
    num_bin_sum <- num_bin_cds + num_bin_fiveUTR + num_bin_threeUTR
    cds_align <- get_align_function(num_bin_cds, trans_info, 
        cds_by_tx0)
    fiveUTR_align <- get_align_function(num_bin_fiveUTR, trans_info, 
        fiveUTR_tx0)
    threeUTR_align <- get_align_function(num_bin_threeUTR, trans_info, 
        threeUTR_tx0)
    align_mtr_5cds3 <- data.frame(fiveUTR_align[, paste0("coordinate.", 
        1:num_bin_fiveUTR)], cds_align[, paste0("coordinate.", 
        1:num_bin_cds)], threeUTR_align[, paste0("coordinate.", 
        1:num_bin_threeUTR)])
    names(align_mtr_5cds3) <- paste0("coordinate.", 1:num_bin_sum)
    width_mtr_5cds3 <- data.frame(coordinate = cbind(replicate(num_bin_fiveUTR, 
        fiveUTR_align[, "width"]/num_bin_fiveUTR), replicate(num_bin_cds, 
        cds_align[, "width"]/num_bin_cds), replicate(num_bin_threeUTR, 
        threeUTR_align[, "width"]/num_bin_threeUTR)))
    if (includeNeighborDNA) {
        num_bin_promoter <- num_bin
        num_bin_tail <- num_bin
        width_promoter <- 1000
        width_tail <- 1000
        num_bin_sum <- num_bin_sum + num_bin_promoter + num_bin_tail
        promoter_align <- get_promoter_function(width_promoter, 
            num_bin_promoter, fiveUTR_align, fiveUTR_tx0, sample_methyl)
        tail_align <- get_tail_function(width_tail, num_bin_tail, 
            threeUTR_align, threeUTR_tx0, sample_methyl)
        align_mtr <- data.frame(coordinate = cbind(promoter_align[, 
            paste0("coordinate.", 1:num_bin_promoter)], fiveUTR_align[, 
            paste0("coordinate.", 1:num_bin_fiveUTR)], cds_align[, 
            paste0("coordinate.", 1:num_bin_cds)], threeUTR_align[, 
            paste0("coordinate.", 1:num_bin_threeUTR)], tail_align[, 
            paste0("coordinate.", 1:num_bin_tail)]))
        names(align_mtr) <- paste0("coordinate.", 1:num_bin_sum)
        width_mtr <- data.frame(coordinate = cbind(replicate(num_bin_promoter, 
            promoter_align[, "width"]/num_bin_promoter), replicate(num_bin_fiveUTR, 
            fiveUTR_align[, "width"]/num_bin_fiveUTR), replicate(num_bin_cds, 
            cds_align[, "width"]/num_bin_cds), replicate(num_bin_threeUTR, 
            threeUTR_align[, "width"]/num_bin_threeUTR), replicate(num_bin_tail, 
            tail_align[, "width"]/num_bin_tail)))
    }
    else {
        align_mtr <- align_mtr_5cds3
        width_mtr <- width_mtr_5cds3
    }
    return(list(align_mtr, width_mtr, trans_info))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
