\name{metaTXplot}
\alias{metaTXplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
metaTXplot(features, txdb, num_bin = 10, includeNeighborDNA = FALSE, comparison = FALSE, lambda = 1, adjust = 0.3, title = "Distribution on mRNA", legend = "sample")
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
  \item{comparison}{
%%     ~~Describe \code{comparison} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{adjust}{
%%     ~~Describe \code{adjust} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
}
  \item{legend}{
%%     ~~Describe \code{legend} here~~
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
function (features, txdb, num_bin = 10, includeNeighborDNA = FALSE, 
    comparison = FALSE, lambda = 1, adjust = 0.3, title = "Distribution on mRNA", 
    legend = "sample") 
{
    print("It may take a few minutes")
    remap_results <- remapCoord(features, txdb, num_bin, includeNeighborDNA)
    align_mtr <- remap_results[[1]]
    width_mtr <- remap_results[[2]]
    trans_info <- remap_results[[3]]
    num_bin_sum <- ncol(align_mtr)
    if (includeNeighborDNA) {
        weight_start <- num_bin_sum/5 + 1
        weight_end <- num_bin_sum * 4/5
        weight_mtr <- replicate(num_bin_sum, rowSums(width_mtr[, 
            weight_start:weight_end]))
        alpha <- get_correct_prob_function(num_bin_sum, align_mtr, 
            weight_mtr^lambda, trans_info)[[1]]
    }
    else {
        weight_mtr <- replicate(num_bin_sum, rowSums(width_mtr))
        alpha <- get_correct_prob_function(num_bin_sum, align_mtr, 
            weight_mtr^lambda, trans_info)[[1]]
    }
    if (!comparison) {
        coord <- 1:num_bin_sum - 0.5
        temp <- data.frame(coord = coord, value = alpha/colSums(width_mtr)/sum(alpha/colSums(width_mtr)), 
            type = legend)
        colnames(temp) <- c("coord", "value", "type")
        data_plot <- temp
        if (includeNeighborDNA) {
            p1 <- ggplot(data_plot, aes(x = coord, group = type, 
                weight = value)) + ggtitle(title) + theme(panel.background = element_blank(), 
                panel.grid.major = element_line(colour = "grey", 
                  linetype = 9, size = 0.2), axis.text.x = element_blank(), 
                axis.ticks = element_blank(), line = element_line(colour = "white", 
                  size = 0.5, linetype = 2, lineend = "butt"), 
                legend.position = "bottom") + geom_density(adjust = adjust, 
                aes(fill = factor(type), colour = factor(type)), 
                alpha = 0.25, fill = "seagreen3", colour = "seagreen3") + 
                xlab("") + ylab("Density\n                     ") + 
                annotate("text", x = num_bin_sum/10, y = -0.007, 
                  label = "Promoter") + annotate("text", x = num_bin_sum * 
                3/10, y = -0.007, label = "5'UTR") + annotate("text", 
                x = num_bin_sum/2, y = -0.007, label = "CDS") + 
                annotate("text", x = num_bin_sum * 7/10, y = -0.007, 
                  label = "3'UTR") + annotate("text", x = num_bin_sum * 
                9/10, y = -0.007, label = "Tail") + geom_vline(xintercept = c(1, 
                2, 3, 4) * num_bin_sum/5, linetype = "dotted", 
                size = 0.5, colour = "black") + annotate("rect", 
                xmin = num_bin_sum * 2/10, xmax = num_bin_sum * 
                  4/10, ymin = -0.0032, ymax = -0.0025, alpha = 0.8, 
                colour = "black") + annotate("rect", xmin = num_bin_sum * 
                4/10, xmax = num_bin_sum * 6/10, ymin = -0.0042, 
                ymax = -0.0017, alpha = 0.3, colour = "black") + 
                annotate("rect", xmin = num_bin_sum * 6/10, xmax = num_bin_sum * 
                  8/10, ymin = -0.0032, ymax = -0.0025, alpha = 0.8, 
                  colour = "black")
        }
    }
    else {
        p1 <- ggplot(data_plot, aes(x = coord, group = type, 
            weight = value)) + ggtitle(title) + theme(panel.background = element_blank(), 
            panel.grid.major = element_line(colour = "grey", 
                linetype = 9, size = 0.2), axis.text.x = element_blank(), 
            axis.ticks = element_blank(), line = element_line(colour = "white", 
                size = 0.5, linetype = 2, lineend = "butt"), 
            legend.position = "bottom") + geom_density(adjust = adjust, 
            aes(fill = factor(type), colour = factor(type)), 
            alpha = 0.25, fill = "seagreen3", colour = "seagreen3") + 
            xlab("") + ylab("Density\n                 ") + annotate("text", 
            x = num_bin_sum * 1/6, y = -0.007, label = "5'UTR") + 
            annotate("text", x = num_bin_sum/2, y = -0.007, label = "CDS") + 
            annotate("text", x = num_bin_sum * 5/6, y = -0.007, 
                label = "3'UTR") + geom_vline(xintercept = c(1, 
            2) * num_bin_sum/3, linetype = "dotted", size = 0.5, 
            colour = "black") + annotate("rect", xmin = num_bin_sum * 
            0, xmax = num_bin_sum * 2/6, ymin = -0.0032, ymax = -0.0025, 
            alpha = 0.8, colour = "black") + annotate("rect", 
            xmin = num_bin_sum * 2/6, xmax = num_bin_sum * 4/6, 
            ymin = -0.0042, ymax = -0.0017, alpha = 0.3, colour = "black") + 
            annotate("rect", xmin = num_bin_sum * 4/6, xmax = num_bin_sum * 
                1, ymin = -0.0032, ymax = -0.0025, alpha = 0.8, 
                colour = "black")
    }
    if (comparison) {
        coord <- 1:num_bin_sum - 0.5
        temp1 <- data.frame(coord = coord, value = colSums(align_mtr)/colSums(width_mtr)/sum(colSums(align_mtr)/colSums(width_mtr)), 
            type = paste0(legend, "_BC"))
        colnames(temp1) <- c("coord", "value", "type")
        temp2 <- data.frame(coord = coord, value = alpha/colSums(width_mtr)/sum(alpha/colSums(width_mtr)), 
            type = paste0(legend, "_AC"))
        colnames(temp2) <- c("coord", "value", "type")
        data_plot <- rbind(temp1, temp2)
        if (includeNeighborDNA) {
            p1 <- ggplot(data_plot, aes(x = coord, group = type, 
                weight = value)) + ggtitle(title) + theme(panel.background = element_blank(), 
                panel.grid.major = element_line(colour = "grey", 
                  linetype = 9, size = 0.2), axis.text.x = element_blank(), 
                axis.ticks = element_blank(), line = element_line(colour = "white", 
                  size = 0.5, linetype = 2, lineend = "butt"), 
                legend.position = "bottom") + geom_density(adjust = adjust, 
                aes(fill = factor(type), colour = factor(type)), 
                alpha = 0.25) + xlab("") + ylab("Density\n                     ") + 
                annotate("text", x = num_bin_sum/10, y = -0.007, 
                  label = "Promoter") + annotate("text", x = num_bin_sum * 
                3/10, y = -0.007, label = "5'UTR") + annotate("text", 
                x = num_bin_sum/2, y = -0.007, label = "CDS") + 
                annotate("text", x = num_bin_sum * 7/10, y = -0.007, 
                  label = "3'UTR") + annotate("text", x = num_bin_sum * 
                9/10, y = -0.007, label = "Tail") + geom_vline(xintercept = c(1, 
                2, 3, 4) * num_bin_sum/5, linetype = "dotted", 
                size = 0.5, colour = "black") + annotate("rect", 
                xmin = num_bin_sum * 2/10, xmax = num_bin_sum * 
                  4/10, ymin = -0.0032, ymax = -0.0025, alpha = 0.8, 
                colour = "black") + annotate("rect", xmin = num_bin_sum * 
                4/10, xmax = num_bin_sum * 6/10, ymin = -0.0042, 
                ymax = -0.0017, alpha = 0.3, colour = "black") + 
                annotate("rect", xmin = num_bin_sum * 6/10, xmax = num_bin_sum * 
                  8/10, ymin = -0.0032, ymax = -0.0025, alpha = 0.8, 
                  colour = "black")
        }
        else {
            p1 <- ggplot(data_plot, aes(x = coord, group = type, 
                weight = value)) + ggtitle(title) + theme(panel.background = element_blank(), 
                panel.grid.major = element_line(colour = "grey", 
                  linetype = 9, size = 0.2), axis.text.x = element_blank(), 
                axis.ticks = element_blank(), line = element_line(colour = "white", 
                  size = 0.5, linetype = 2, lineend = "butt"), 
                legend.position = "bottom") + geom_density(adjust = adjust, 
                aes(fill = factor(type), colour = factor(type)), 
                alpha = 0.25) + xlab("") + ylab("Density\n                     ") + 
                annotate("text", x = num_bin_sum * 1/6, y = -0.007, 
                  label = "5'UTR") + annotate("text", x = num_bin_sum/2, 
                y = -0.007, label = "CDS") + annotate("text", 
                x = num_bin_sum * 5/6, y = -0.007, label = "3'UTR") + 
                geom_vline(xintercept = c(1, 2) * num_bin_sum/3, 
                  linetype = "dotted", size = 0.5, colour = "black") + 
                annotate("rect", xmin = num_bin_sum * 0, xmax = num_bin_sum * 
                  2/6, ymin = -0.0032, ymax = -0.0025, alpha = 0.8, 
                  colour = "black") + annotate("rect", xmin = num_bin_sum * 
                2/6, xmax = num_bin_sum * 4/6, ymin = -0.0042, 
                ymax = -0.0017, alpha = 0.3, colour = "black") + 
                annotate("rect", xmin = num_bin_sum * 4/6, xmax = num_bin_sum * 
                  1, ymin = -0.0032, ymax = -0.0025, alpha = 0.8, 
                  colour = "black")
        }
        return(p1)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
