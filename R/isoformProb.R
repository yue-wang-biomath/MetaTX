isoformProb <-
function(remap_results,
        includeNeighborDNA = FALSE,
        lambda = 2     
    ){
        align_mtr          <- remap_results[[1]]
        width_mtr          <- remap_results[[2]]
        trans_info         <- remap_results[[3]]
        num_bin_sum        <- ncol(align_mtr)
        
        get_correct_prob_function <-
          function(num_bin_sum, align_mtr, weight_mtr, trans_info){
            
            alpha                     <- matrix(1/num_bin_sum, 1, num_bin_sum)
            index_methyl              <- trans_info[, 'index_methyl']
            
            for (j in 1:20){
              numerator_prob          <- align_mtr 
              
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
        
        if(includeNeighborDNA){
            weight_mtr   <- replicate(num_bin_sum, rowSums(width_mtr))
            alpha        <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[1]]
            prob_mtr     <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[2]]
        }else{
            weight_start <- num_bin_sum / 5 + 1
            weight_end   <- num_bin_sum * 4/ 5
            weight_mtr   <- replicate(num_bin_sum, rowSums(width_mtr[, weight_start:weight_end]))    
            alpha        <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[1]]
            prob_mtr     <- get_correct_prob_function(num_bin_sum, align_mtr, weight_mtr^lambda, trans_info)[[2]]
        } 
        prob           <- as.numeric(rowSums(prob_mtr))
        prob_isoforms  <-  cbind(trans_info,
                                 isoform_prob = prob)
        return(prob_isoforms)
    }
