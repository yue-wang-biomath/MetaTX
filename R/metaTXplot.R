metaTXplot <-
function(    remap_results,
             num_bin = 10,
             includeNeighborDNA = FALSE,
             comparison = FALSE,
             lambda = 2,
             adjust = 0.15,
             title  = 'Distribution on mRNA',
             legend = 'sample'){

    # function 1
    get_correct_prob_function <-
    function(num_bin_sum, align_mtr, weight_mtr, trans_info, lambda){
      
      alpha                     <- matrix(1/num_bin_sum, 1, num_bin_sum)
      index_methyl              <- trans_info[, 'index_methyl']
      
      for (j in 1:20){
        numerator_prob          <- align_mtr * width_mtr
        
        for (k in 1:num_bin_sum){
          numerator_prob[, k]   <- numerator_prob[, k] * alpha[k]
        }
        row_sum_prob            <- as.vector(rowSums(numerator_prob))
        row_sum_prob            <- data.frame(group_name = as.character(index_methyl), value = row_sum_prob)
        sum_prob                <- aggregate(row_sum_prob[,'value'], 
                                             by = list(group_name = factor(trans_info[,'index_methyl'], levels = unique(trans_info[,'index_methyl']))), 
                                             FUN = sum)
        names(sum_prob)         <- c('index_methyl', 'value')
        denominator_prob        <- sum_prob[, 'value']
        names(denominator_prob) <- sum_prob[, 'index_methyl']
        denominator_prob        <- denominator_prob[as.character(trans_info[, 'index_methyl'])]
        denominator_prob        <- replicate(num_bin_sum, denominator_prob)
        denominator_prob[denominator_prob == 0] <- 1
        
        prob_mtr                <- numerator_prob / denominator_prob 
        alpha_numerator         <- colSums(prob_mtr * weight_mtr)
        alpha                   <- alpha_numerator/sum(alpha_numerator)
      }
      
      
      return(list(alpha, prob_mtr))
    }

        align_mtr          <- remap_results[[1]]
        width_mtr          <- remap_results[[2]]
        trans_info         <- remap_results[[3]]
        num_bin_sum        <- ncol(align_mtr)
        
        if(includeNeighborDNA){
            weight_start     <- num_bin_sum / 5 + 1
            weight_end       <- num_bin_sum * 4/ 5
            weight_mtr       <- replicate(num_bin_sum, rowSums(width_mtr[, weight_start:weight_end]))    
            alpha            <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[1]]
        }else{
            weight_mtr    <- replicate(num_bin_sum, rowSums(width_mtr))
            alpha         <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[1]]
        } 
        
        
        
        if(!comparison){
            coord        <- 1:num_bin_sum - 0.5
            temp         <- data.frame(coord = coord
                                       , value =  alpha / colSums(width_mtr) / sum(alpha / colSums(width_mtr))
                                       , type  =  legend)
            colnames(temp)<- c('coord', 'value', 'type')
            data_plot     <- temp
            
            
            if(includeNeighborDNA){
                p1 <- 
                    ggplot(data_plot, aes(x=coord, group=type, weight=value)) +
                    ggtitle(title) +
                    theme(panel.background =element_blank(),
                          panel.grid.major = element_line(colour = 'grey', linetype = 9, size = 0.2),
                          axis.text.x = element_blank(), axis.ticks = element_blank(),
                          line = element_line(colour = "white", size = 0.5, linetype = 2, lineend = "butt"),
                          legend.position = 'bottom') + 
                    geom_density(adjust = adjust, aes(fill=factor(type), colour = factor(type)),alpha = 0.25, fill = 'seagreen3', colour = 'seagreen3') +
                    xlab("") + 
                    ylab("Density
                         ") +
                    annotate("text", x = num_bin_sum / 10, y = -0.007, label = "Promoter") +
                    annotate("text", x = num_bin_sum * 3 / 10, y = -0.007, label = "5'UTR") +
                    annotate("text", x = num_bin_sum / 2, y = -0.007, label = "CDS") +
                    annotate("text", x = num_bin_sum * 7 / 10, y = -0.007, label = "3'UTR") +
                    annotate("text", x = num_bin_sum * 9 / 10, y = -0.007, label = "Tail") + 
                    geom_vline(xintercept= c(1, 2, 3, 4) * num_bin_sum / 5, linetype = "dotted", size = 0.5,colour = 'black') +
                    annotate("rect", xmin = num_bin_sum * 2 / 10, xmax = num_bin_sum * 4 / 10, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")+
                    annotate("rect", xmin = num_bin_sum * 4 / 10, xmax = num_bin_sum * 6 / 10, ymin = -0.0042, ymax = -0.0017, alpha = .3, colour = "black")+
                    annotate("rect", xmin = num_bin_sum * 6 / 10, xmax = num_bin_sum * 8 / 10, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")
            }  
        else{
            p1 <- 
                ggplot(data_plot, aes(x=coord, group=type, weight=value)) +
                ggtitle(title) +
                theme(panel.background =element_blank(),
                      panel.grid.major = element_line(colour = 'grey', linetype = 9, size = 0.2),
                      axis.text.x = element_blank(), axis.ticks = element_blank(),
                      line = element_line(colour = "white", size = 0.5, linetype = 2, lineend = "butt"),
                      legend.position = 'bottom') + 
                geom_density(adjust = adjust, aes(fill=factor(type), colour = factor(type)),alpha = 0.25, fill = 'seagreen3', colour = 'seagreen3') +
                xlab("") + 
                ylab("Density
                     ") +
                annotate("text", x = num_bin_sum * 1/6, y = -0.007, label = "5'UTR") +
                annotate("text", x = num_bin_sum / 2, y = -0.007, label = "CDS") +
                annotate("text", x = num_bin_sum * 5 / 6, y = -0.007, label = "3'UTR") +
                geom_vline(xintercept= c(1, 2) * num_bin_sum /3 , linetype = "dotted", size = 0.5,colour = 'black') +
                annotate("rect", xmin = num_bin_sum * 0,  xmax = num_bin_sum *2 / 6, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")+
                annotate("rect", xmin = num_bin_sum * 2 / 6, xmax = num_bin_sum * 4 / 6, ymin = -0.0042, ymax = -0.0017, alpha = .3, colour = "black")+
                annotate("rect", xmin = num_bin_sum * 4 / 6, xmax = num_bin_sum * 1, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")
        }
       }
        
        if(comparison){
            coord            <- 1:num_bin_sum - 0.5
            temp1            <- data.frame(coord = coord
                                        , value =  colSums(align_mtr)  / colSums(width_mtr)  / sum( colSums(align_mtr)  / colSums(width_mtr) )
                                        , type  = paste0(legend, '_BC'))
            colnames(temp1)  <- c('coord', 'value', 'type')
            row.names(temp1) <- 1:num_bin_sum
            temp2            <- data.frame(coord  = coord
                                        , value =  alpha / colSums(width_mtr) / sum(alpha / colSums(width_mtr))
                                        , type  = paste0(legend, '_AC'))
            colnames(temp2)  <- c('coord', 'value', 'type')
            row.names(temp2) <- 1:num_bin_sum + num_bin_sum
            data_plot        <- rbind(temp1, temp2)
            
            
            if(includeNeighborDNA){
                p1 <- 
                    ggplot(data_plot, aes(x=coord, group=type, weight=value)) +
                    ggtitle(title) +
                    theme(panel.background =element_blank(),
                          panel.grid.major = element_line(colour = 'grey', linetype = 9, size = 0.2),
                          axis.text.x = element_blank(), axis.ticks = element_blank(),
                          line = element_line(colour = "white", size = 0.5, linetype = 2, lineend = "butt"),
                          legend.position = 'bottom') + 
                    geom_density(adjust = adjust, aes(fill=factor(type), colour = factor(type)),alpha = 0.25) +
                    xlab("") + 
                    ylab("Density
                         ") +
                    annotate("text", x = num_bin_sum / 10, y = -0.007, label = "Promoter") +
                    annotate("text", x = num_bin_sum * 3 / 10, y = -0.007, label = "5'UTR") +
                    annotate("text", x = num_bin_sum / 2, y = -0.007, label = "CDS") +
                    annotate("text", x = num_bin_sum * 7 / 10, y = -0.007, label = "3'UTR") +
                    annotate("text", x = num_bin_sum * 9 / 10, y = -0.007, label = "Tail") + 
                    geom_vline(xintercept= c(1, 2, 3, 4) * num_bin_sum / 5, linetype = "dotted", size = 0.5,colour = 'black') +
                    annotate("rect", xmin = num_bin_sum * 2 / 10, xmax = num_bin_sum * 4 / 10, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")+
                    annotate("rect", xmin = num_bin_sum * 4 / 10, xmax = num_bin_sum * 6 / 10, ymin = -0.0042, ymax = -0.0017, alpha = .3, colour = "black")+
                    annotate("rect", xmin = num_bin_sum * 6 / 10, xmax = num_bin_sum * 8 / 10, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")
                
            }else{
                p1 <- 
                    ggplot(data_plot, aes(x=coord, group=type, weight=value)) +
                    ggtitle(title) +
                    theme(panel.background = element_blank(),
                          panel.grid.major = element_line(colour = 'grey', linetype = 9, size = 0.2),
                          axis.text.x = element_blank(), axis.ticks = element_blank(),
                          line = element_line(colour = "white", size = 0.5, linetype = 2, lineend = "butt"),
                          legend.position = 'bottom') + 
                    geom_density(adjust = adjust, aes(fill=factor(type), colour = factor(type)),alpha = 0.25) +
                    xlab("") + 
                    ylab("Density
                         ") +
                    annotate("text", x = num_bin_sum * 1/ 6, y = -0.007, label = "5'UTR") +
                    annotate("text", x = num_bin_sum / 2, y = -0.007, label = "CDS") +
                    annotate("text", x = num_bin_sum * 5 / 6, y = -0.007, label = "3'UTR") +
                    geom_vline(xintercept= c(1, 2) * num_bin_sum / 3 , linetype = "dotted", size = 0.5,colour = 'black') +
                    annotate("rect", xmin = num_bin_sum * 0,  xmax = num_bin_sum * 2 / 6, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")+
                    annotate("rect", xmin = num_bin_sum * 2 / 6, xmax = num_bin_sum * 4 / 6, ymin = -0.0042, ymax = -0.0017, alpha = .3, colour = "black")+
                    annotate("rect", xmin = num_bin_sum * 4 / 6, xmax = num_bin_sum * 1, ymin = -0.0032, ymax = -0.0025, alpha = .8, colour = "black")
            }
            
        }
        return(p1)  
    }
