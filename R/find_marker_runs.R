find_marker_runs <- function(
    markers,
    VECTOR_ALLOC_CHUNK_SIZE = 200000,
    MIN_RUN_LENGTH = 4,
    MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE = 0.7,
    CNV_MIN_DELTA_ABS = 0.2,
    DIPLOID_MIN_DELTA_ABS = 0.3
    ){
        

    is_potential_run_active <- FALSE
    run_chrom <- vector("character", VECTOR_ALLOC_CHUNK_SIZE)
    run_start <- vector("integer", VECTOR_ALLOC_CHUNK_SIZE)
    run_end <- vector("integer", VECTOR_ALLOC_CHUNK_SIZE)
    run_i <- 0

    for(i in seq_len(nrow(markers))) {
        chrom_i <- markers$chrom[i]
        pos_i <- markers$pos[i]
        reads.RNA.ref_i <- markers$reads.RNA.ref[i]
        reads.RNA.alt_i <- markers$reads.RNA.alt[i]
        p_value_i <- markers$p_value[i]
        delta_abs_i <- markers$delta_abs[i]
        is_cnv_i <- (markers$copynumber_tag[i] == "cnv")
        
        
        
        is_marker_valid <- (reads.RNA.ref_i + reads.RNA.alt_i) > 5
        if(!is_marker_valid){next}
        
        is_marker_extreme <- ((reads.RNA.ref_i == 0) | (reads.RNA.alt_i == 0) | (((reads.RNA.ref_i <= 1) | (reads.RNA.alt_i <= 1)) & (p_value_i <= 0.05)))
        is_marker_extreme_with_error <- (((reads.RNA.ref_i <= 1) | (reads.RNA.alt_i <= 1)) & (p_value_i > 0.05))
        is_marker_significant <- (p_value_i <= 0.05) & (delta_abs_i >= ifelse(is_cnv_i, CNV_MIN_DELTA_ABS, DIPLOID_MIN_DELTA_ABS))
        is_marker_negative <- (!is_marker_extreme) & (!is_marker_extreme_with_error)
        
        # check for next chromosome or two high distance
        if(is_potential_run_active){
            if(!(temp_chrom == chrom_i)){
                delete_last_marker <- FALSE
                delete_last_two_markers <- FALSE
                
                if(temp_negative[temp_i]){
                    end_run <- TRUE
                    delete_last_marker <- FALSE
                }
                if(temp_i >= 2){
                    if(temp_negative[temp_i - 1]){
                    end_run <- TRUE
                    delete_last_two_markers <- FALSE
                    }
                }
                
                
                if(delete_last_two_markers){
                    #evaluate run
                    if(
                    ((temp_i-2) >= MIN_RUN_LENGTH) &
                    ((sum(temp_extreme[1:(temp_i-2)] | temp_significant[1:(temp_i-2)]) / (temp_i-2)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                    run_i <- run_i +1
                    run_chrom[run_i] <- temp_chrom
                    run_start[run_i] <- temp_pos[1]
                    run_end[run_i] <- temp_pos[temp_i-2]
                    is_potential_run_active <- FALSE
                    }else{
                    # discard run
                    is_potential_run_active <- FALSE
                    }
                }else if (delete_last_marker){
                    #evaluate run
                    if(
                    ((temp_i-1) >= MIN_RUN_LENGTH) &
                    ((sum(temp_extreme[1:(temp_i-1)] | temp_significant[1:(temp_i-1)]) / (temp_i-1))>= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                    run_i <- run_i +1
                    run_chrom[run_i] <- temp_chrom
                    run_start[run_i] <- temp_pos[1]
                    run_end[run_i] <- temp_pos[temp_i-1]
                    is_potential_run_active <- FALSE
                    }else{
                    # discard run
                    is_potential_run_active <- FALSE
                    }
                } else {
                    #evaluate run
                    if(
                    ((temp_i) >= MIN_RUN_LENGTH) &
                    ((sum(temp_extreme[1:(temp_i)] | temp_significant[1:(temp_i)]) / (temp_i)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                    run_i <- run_i +1
                    run_chrom[run_i] <- temp_chrom
                    run_start[run_i] <- temp_pos[1]
                    run_end[run_i] <- temp_pos[temp_i]
                    is_potential_run_active <- FALSE
                    }else{
                    # discard run
                    is_potential_run_active <- FALSE
                    }
                }
            }
        }
        
        if(is_marker_negative){
            if(is_potential_run_active){
                end_run <- FALSE
                delete_last_marker <- FALSE
                delete_last_two_markers <- FALSE
                
                if(temp_negative[temp_i]){
                    end_run <- TRUE
                    delete_last_marker <- FALSE
                }
                if(temp_i >= 2){
                    if(temp_negative[temp_i - 1]){
                    end_run <- TRUE
                    delete_last_two_markers <- FALSE
                    }
                }
                
                if(end_run){
                    if(delete_last_two_markers){
                    #evaluate run
                    if(
                        ((temp_i-2) >= MIN_RUN_LENGTH) &
                        ((sum(temp_extreme[1:(temp_i-2)] | temp_significant[1:(temp_i-2)]) / (temp_i-2)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                        run_i <- run_i +1
                        run_chrom[run_i] <- temp_chrom
                        run_start[run_i] <- temp_pos[1]
                        run_end[run_i] <- temp_pos[temp_i-2]
                        is_potential_run_active <- FALSE
                    }else{
                        # discard run
                        is_potential_run_active <- FALSE
                    }
                    }else if (delete_last_marker){
                    #evaluate run
                    if(
                        ((temp_i-1) >= MIN_RUN_LENGTH) &
                        ((sum(temp_extreme[1:(temp_i-1)] | temp_significant[1:(temp_i-1)]) / (temp_i-1)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                        run_i <- run_i +1
                        run_chrom[run_i] <- temp_chrom
                        run_start[run_i] <- temp_pos[1]
                        run_end[run_i] <- temp_pos[temp_i-1]
                        is_potential_run_active <- FALSE
                    }else{
                        # discard run
                        is_potential_run_active <- FALSE
                    }
                    } else {
                    # should actually never be the case
                    #evaluate run
                    if(
                        ((temp_i) >= MIN_RUN_LENGTH) &
                        ((sum(temp_extreme[1:(temp_i)] | temp_significant[1:(temp_i)]) / (temp_i)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                        run_i <- run_i +1
                        run_chrom[run_i] <- temp_chrom
                        run_start[run_i] <- temp_pos[1]
                        run_end[run_i] <- temp_pos[temp_i]
                        is_potential_run_active <- FALSE
                    }else{
                        # discard run
                        is_potential_run_active <- FALSE
                    }
                    }
                }else{
                    # add to run
                    temp_i <- temp_i +1
                    temp_pos[temp_i] <- pos_i
                    temp_extreme[temp_i] <- is_marker_extreme
                    temp_extreme_with_error[temp_i] <- is_marker_extreme_with_error
                    temp_significant[temp_i] <- is_marker_significant
                    temp_negative[temp_i] <- is_marker_negative
                }
            } else {
            # do nothing
            1
            }
        } else {
            if(is_potential_run_active){
                # add to run
                temp_i <- temp_i +1
                temp_pos[temp_i] <- pos_i
                temp_extreme[temp_i] <- is_marker_extreme
                temp_extreme_with_error[temp_i] <- is_marker_extreme_with_error
                temp_significant[temp_i] <- is_marker_significant
                temp_negative[temp_i] <- is_marker_negative
            }else{
                # start run
                temp_chrom <- chrom_i
                temp_pos <- vector("integer", VECTOR_ALLOC_CHUNK_SIZE)
                temp_extreme <- vector("logical", VECTOR_ALLOC_CHUNK_SIZE)
                temp_extreme_with_error <- vector("logical", VECTOR_ALLOC_CHUNK_SIZE)
                temp_significant <- vector("logical", VECTOR_ALLOC_CHUNK_SIZE)
                temp_negative <- vector("logical", VECTOR_ALLOC_CHUNK_SIZE)
                temp_i <- 0
                is_potential_run_active <- TRUE
                
                # add to run
                temp_i <- temp_i +1
                temp_pos[temp_i] <- pos_i
                temp_extreme[temp_i] <- is_marker_extreme
                temp_extreme_with_error[temp_i] <- is_marker_extreme_with_error
                temp_significant[temp_i] <- is_marker_significant
                temp_negative[temp_i] <- is_marker_negative
            }
            
        }
        
        if((i == nrow(markers)) & is_potential_run_active){
            # evaluate run
            delete_last_marker <- FALSE
            delete_last_two_markers <- FALSE
            
            if(temp_negative[temp_i]){
                end_run <- TRUE
                delete_last_marker <- FALSE
                }
                if(temp_i >= 2){
                    if(temp_negative[temp_i - 1]){
                        end_run <- TRUE
                        delete_last_two_markers <- FALSE
                    }
            }
            
            
            if(delete_last_two_markers){
                #evaluate run
                if(
                    ((temp_i-2) >= MIN_RUN_LENGTH) &
                    ((sum(temp_extreme[1:(temp_i-2)] | temp_significant[1:(temp_i-2)]) / (temp_i-2)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                ){
                    run_i <- run_i +1
                    run_chrom[run_i] <- temp_chrom
                    run_start[run_i] <- temp_pos[1]
                    run_end[run_i] <- temp_pos[temp_i-2]
                    is_potential_run_active <- FALSE
                }else{
                    # discard run
                    is_potential_run_active <- FALSE
                }
                }else if (delete_last_marker){
                    #evaluate run
                    if(
                        ((temp_i-1) >= MIN_RUN_LENGTH) &
                        ((sum(temp_extreme[1:(temp_i-1)] | temp_significant[1:(temp_i-1)]) / (temp_i-1)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                        run_i <- run_i +1
                        run_chrom[run_i] <- temp_chrom
                        run_start[run_i] <- temp_pos[1]
                        run_end[run_i] <- temp_pos[temp_i-1]
                        is_potential_run_active <- FALSE
                    }else{
                        # discard run
                        is_potential_run_active <- FALSE
                    }
                } else {
                #evaluate run
                    if(
                        ((temp_i) >= MIN_RUN_LENGTH) &
                        ((sum(temp_extreme[1:(temp_i)] | temp_significant[1:(temp_i)]) / (temp_i)) >= MIN_EXTREME_OR_SIGNIFIGANT_PERCENTAGE)
                    ){
                        run_i <- run_i +1
                        run_chrom[run_i] <- temp_chrom
                        run_start[run_i] <- temp_pos[1]
                        run_end[run_i] <- temp_pos[temp_i]
                        is_potential_run_active <- FALSE
                    }else{
                        # discard run
                        is_potential_run_active <- FALSE
                    }
            }
        }
    }

    runs <- data.table::setDT(data.frame(
        chrom = run_chrom[1:run_i],
        start = run_start[1:run_i],
        end = run_end[1:run_i]
    ))
    return(runs)
}