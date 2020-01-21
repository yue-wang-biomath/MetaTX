metaTXcount <-
function(features,
             txdb,
             num_bin = 10,
             includeNeighborDNA = FALSE
    ){
        remap_results      <- remapCoord(methyl, txdb, num_bin, includeNeighborDNA)
        align_mtr          <- remap_results[[1]]
        num_bin_sum        <- ncol(align_mtr)
        feature_counts     <- colSums(align_mtr)
        return(feature_counts)
    }
