metaTXcount <-
function(remap_results
    ){
        align_mtr          <- remap_results[[1]]
        num_bin_sum        <- ncol(align_mtr)
        feature_counts     <- colSums(align_mtr)
        return(feature_counts)
    }
