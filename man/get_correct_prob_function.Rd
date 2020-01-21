\name{get_correct_prob_function}
\alias{get_correct_prob_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample data for checking CDS region
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr, trans_info)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{num_bin_sum}{
%%     ~~Describe \code{num_bin_sum} here~~
}
  \item{align_mtr}{
%%     ~~Describe \code{align_mtr} here~~
}
  \item{weight_mtr}{
%%     ~~Describe \code{weight_mtr} here~~
}
  \item{trans_info}{
%%     ~~Describe \code{trans_info} here~~
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
function (num_bin_sum, align_mtr, weight_mtr, trans_info) 
{
    alpha <- matrix(1/num_bin_sum, 1, num_bin_sum)
    index_methyl <- trans_info[, "index_methyl"]
    for (j in 1:20) {
        numerator_prob <- align_mtr
        for (k in 1:num_bin_sum) {
            numerator_prob[, k] <- numerator_prob[, k] * alpha[k]
        }
        row_sum_prob <- as.vector(rowSums(numerator_prob))
        row_sum_prob <- data.frame(group_name = as.character(index_methyl), 
            value = row_sum_prob)
        sum_prob <- aggregate(row_sum_prob[, "value"], by = list(group_name = factor(trans_info[, 
            "index_methyl"], levels = unique(trans_info[, "index_methyl"]))), 
            FUN = sum)
        names(sum_prob) <- c("index_methyl", "value")
        denominator_prob <- sum_prob[, "value"]
        names(denominator_prob) <- sum_prob[, "index_methyl"]
        denominator_prob <- denominator_prob[as.character(trans_info[, 
            "index_methyl"])]
        denominator_prob <- replicate(num_bin_sum, denominator_prob)
        denominator_prob[denominator_prob == 0] <- 1
        prob_mtr <- numerator_prob/denominator_prob
        alpha_numerator <- colSums(prob_mtr * weight_mtr)
        alpha <- alpha_numerator/sum(alpha_numerator)
    }
    return(list(alpha, prob_mtr))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
