evaluate_variant_position_by_single_location <- function(
    variant_chrom,
    variant_position,
    gene_chrom,
    gene_start,
    gene_end,
    reg_region_chrom,
    reg_region_start,
    reg_region_end,
    gene_TAD_start,
    gene_TAD_end,
    reg_region_TAD_start,
    reg_region_TAD_end
    ) {
    same_chrom_gene_and_reg_region <- (gene_chrom == reg_region_chrom)
    same_chrom_gene_and_variant <- (gene_chrom == variant_chrom)
    same_chrom_variant_and_reg_region <- (variant_chrom == reg_region_chrom)

    same_TAD_gene_and_reg_region <- same_chrom_gene_and_reg_region &
        (reg_region_end >= gene_TAD_start) &
        (reg_region_start <= gene_TAD_end)
    same_TAD_gene_and_variant <- same_chrom_gene_and_variant &
        (   #within TAD
            (variant_position >= gene_TAD_start) &
            (variant_position <= gene_TAD_end)
        ) |
        ( # OR within gene
            (variant_position >= gene_start) &
            (variant_position <= gene_end)
        )
    same_TAD_variant_and_reg_region <- same_chrom_variant_and_reg_region &
        (   #within TAD
            (variant_position >= reg_region_TAD_start) &
            (variant_position <= reg_region_TAD_end)
        ) |
        ( # within genehancer
            (variant_position >= reg_region_start) &
            (variant_position <= reg_region_end)
        )

    gene_and_reg_region_overlap <- same_chrom_gene_and_reg_region &
        (reg_region_end >= gene_start) &
        (reg_region_start <= gene_end)

    gene_and_reg_region_adjacent <- same_chrom_gene_and_reg_region &
        (
            ((reg_region_end + 1) == gene_start) |
            ((gene_end +1 ) == reg_region_start)
        )

    start_region_between_gene_and_reg_region <- (
        ifelse(
            (!same_chrom_gene_and_reg_region) |
            gene_and_reg_region_overlap |
            gene_and_reg_region_adjacent
            ,
            NA,
            ifelse(
                # reg region positioned before gene
                ((reg_region_end + 1) <= (gene_start - 1)),
                reg_region_end + 1,
                ifelse(
                    # gene positioned before reg region
                    ((gene_end + 1) <= (reg_region_start - 1)),
                    gene_end + 1,
                    NA
                )
            )
        )
    )

    end_region_between_gene_and_reg_region <- (
        ifelse(
            (!same_chrom_gene_and_reg_region) |
            gene_and_reg_region_overlap |
            gene_and_reg_region_adjacent
            ,
            NA,
            ifelse(
                # reg region positioned before gene
                ((reg_region_end + 1) <= (gene_start - 1)),
                gene_start - 1,
                ifelse(
                    # gene positioned before reg region
                    ((gene_end + 1) <= (reg_region_start - 1)),
                    reg_region_start - 1,
                    NA
                )
            )
        )
    )

    position_category <- ifelse(
        same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
        same_TAD_gene_and_variant & same_TAD_variant_and_reg_region & # in shared TAD
        (variant_position > gene_end) & (variant_position > reg_region_end), # high of both
        "ON_SHARED_CHROM__IN_SHARED_TAD__HIGH_OF_BOTH",
        ifelse(
            same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
            same_TAD_gene_and_variant & same_TAD_variant_and_reg_region & # in shared TAD
            (variant_position < gene_start) & (variant_position < reg_region_start), # low of both
            "ON_SHARED_CHROM__IN_SHARED_TAD__LOW_OF_BOTH",
            ifelse(
                same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                same_TAD_gene_and_variant & same_TAD_variant_and_reg_region & # in shared TAD
                (!is.na(start_region_between_gene_and_reg_region)) & 
                (!is.na(end_region_between_gene_and_reg_region)) &
                (!gene_and_reg_region_overlap) &
                (!gene_and_reg_region_adjacent) & # between region exists (duplicate conditions to catch false data)
                (variant_position >= start_region_between_gene_and_reg_region) &
                (variant_position <= end_region_between_gene_and_reg_region), # within between region
                "ON_SHARED_CHROM__IN_SHARED_TAD__BETWEEN",
                ifelse(
                    same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                    same_TAD_gene_and_variant & same_TAD_variant_and_reg_region & # in shared TAD
                    (variant_position >= reg_region_start) & (variant_position <= reg_region_end) & # within genehancer
                    (!((variant_position >= gene_start) & (variant_position <= gene_end))), # not within gene
                    "ON_SHARED_CHROM__IN_SHARED_TAD__IN_REG",
                    ifelse(
                        same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                        same_TAD_gene_and_variant & same_TAD_variant_and_reg_region & # in shared TAD
                        (variant_position >= gene_start) & (variant_position <= gene_end) & # within gene
                        (!((variant_position >= reg_region_start) & (variant_position <= reg_region_end))), # not within genehancer
                        "ON_SHARED_CHROM__IN_SHARED_TAD__IN_GENE",
                        ifelse(
                            same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                            same_TAD_gene_and_variant & same_TAD_variant_and_reg_region & # in shared TAD
                            gene_and_reg_region_overlap &
                            (variant_position >= gene_start) & (variant_position <= gene_end) & # within gene
                            (variant_position >= reg_region_start) & (variant_position <= reg_region_end), #within genehancer
                            "ON_SHARED_CHROM__IN_SHARED_TAD__IN_OVERLAP",
                            ifelse(
                                same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                (!same_TAD_gene_and_variant) & same_TAD_variant_and_reg_region & # in reg region TAD
                                (variant_position > gene_end) & (variant_position > reg_region_end), # high of both
                                "ON_SHARED_CHROM__IN_REG_TAD__HIGH_OF_BOTH",
                                ifelse(
                                    same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                    (!same_TAD_gene_and_variant) & same_TAD_variant_and_reg_region & # in reg region TAD
                                    (variant_position < gene_start) & (variant_position < reg_region_start), # low of both
                                    "ON_SHARED_CHROM__IN_REG_TAD__LOW_OF_BOTH",
                                    ifelse(
                                        same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                        (!same_TAD_gene_and_variant) & same_TAD_variant_and_reg_region & # in reg region TAD
                                        (!is.na(start_region_between_gene_and_reg_region)) & 
                                        (!is.na(end_region_between_gene_and_reg_region)) &
                                        (!gene_and_reg_region_overlap) &
                                        (!gene_and_reg_region_adjacent) & # between region exists (duplicate conditions to catch false data)
                                        (variant_position >= start_region_between_gene_and_reg_region) &
                                        (variant_position <= end_region_between_gene_and_reg_region), # within between region
                                        "ON_SHARED_CHROM__IN_REG_TAD__BETWEEN",
                                        ifelse(
                                            same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                            (!same_TAD_gene_and_variant) & same_TAD_variant_and_reg_region & # in reg region TAD
                                            (variant_position >= reg_region_start) & (variant_position <= reg_region_end) & # within genehancer
                                            (!((variant_position >= gene_start) & (variant_position <= gene_end))), # not within gene
                                            "ON_SHARED_CHROM__IN_REG_TAD__IN_REG",
                                            ifelse(
                                                same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                                same_TAD_gene_and_variant & (!same_TAD_variant_and_reg_region) & # in gene TAD
                                                (variant_position > gene_end) & (variant_position > reg_region_end), # high of both
                                                "ON_SHARED_CHROM__IN_GENE_TAD__HIGH_OF_BOTH",
                                                ifelse(
                                                    same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                                    same_TAD_gene_and_variant & (!same_TAD_variant_and_reg_region) & # in gene TAD
                                                    (variant_position < gene_start) & (variant_position < reg_region_start), # low of both
                                                    "ON_SHARED_CHROM__IN_GENE_TAD__LOW_OF_BOTH",
                                                    ifelse(
                                                        same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                                        same_TAD_gene_and_variant & (!same_TAD_variant_and_reg_region) & # in gene TAD
                                                        (!is.na(start_region_between_gene_and_reg_region)) & 
                                                        (!is.na(end_region_between_gene_and_reg_region)) &
                                                        (!gene_and_reg_region_overlap) &
                                                        (!gene_and_reg_region_adjacent) & # between region exists (duplicate conditions to catch false data)
                                                        (variant_position >= start_region_between_gene_and_reg_region) &
                                                        (variant_position <= end_region_between_gene_and_reg_region), # within between region
                                                        "ON_SHARED_CHROM__IN_GENE_TAD__BETWEEN",
                                                        ifelse(
                                                            same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                                            same_TAD_gene_and_variant & (!same_TAD_variant_and_reg_region) & # in gene TAD
                                                            (variant_position >= gene_start) & (variant_position <= gene_end) & # within gene
                                                            (!((variant_position >= reg_region_start) & (variant_position <= reg_region_end))), # not within reg region
                                                            "ON_SHARED_CHROM__IN_GENE_TAD__IN_GENE",
                                                            ifelse(
                                                                same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                                                (!same_TAD_gene_and_variant) & (!same_TAD_variant_and_reg_region) & # outside TAD
                                                                (variant_position > gene_end) & (variant_position > reg_region_end), # high of both
                                                                "ON_SHARED_CHROM__OUTSIDE_TAD__HIGH_OF_BOTH",
                                                                ifelse(
                                                                    same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                                                    (!same_TAD_gene_and_variant) & (!same_TAD_variant_and_reg_region) & # outside TAD
                                                                    (variant_position < gene_start) & (variant_position < reg_region_start), # low of both
                                                                    "ON_SHARED_CHROM__OUTSIDE_TAD__LOW_OF_BOTH",
                                                                    ifelse(
                                                                        same_chrom_gene_and_reg_region & same_chrom_gene_and_variant & same_chrom_variant_and_reg_region & # on shared chrom
                                                                        (!same_TAD_gene_and_variant) & (!same_TAD_variant_and_reg_region) & # outside TAD
                                                                        (!is.na(start_region_between_gene_and_reg_region)) & 
                                                                        (!is.na(end_region_between_gene_and_reg_region)) &
                                                                        (!gene_and_reg_region_overlap) &
                                                                        (!gene_and_reg_region_adjacent) & # between region exists (duplicate conditions to catch false data)
                                                                        (variant_position >= start_region_between_gene_and_reg_region) &
                                                                        (variant_position <= end_region_between_gene_and_reg_region), # within between region
                                                                        "ON_SHARED_CHROM__OUTSIDE_TAD__BETWEEN",
                                                                        ifelse(
                                                                            (!same_chrom_gene_and_reg_region) & (!same_chrom_gene_and_variant) & same_chrom_variant_and_reg_region & # on reg region chrom
                                                                            same_TAD_variant_and_reg_region & # in reg region TAD
                                                                            (variant_position >= reg_region_start) & (variant_position <= reg_region_end), # within reg region
                                                                            "ON_REG_CHROM__IN_REG_TAD__IN_REG",
                                                                            ifelse(
                                                                                (!same_chrom_gene_and_reg_region) & (!same_chrom_gene_and_variant) & same_chrom_variant_and_reg_region & # on reg region chrom
                                                                                same_TAD_variant_and_reg_region & # in reg region TAD
                                                                                (variant_position > reg_region_end), # high of reg region
                                                                                "ON_REG_CHROM__IN_REG_TAD__HIGH_OF_REG",
                                                                                ifelse(
                                                                                    (!same_chrom_gene_and_reg_region) & (!same_chrom_gene_and_variant) & same_chrom_variant_and_reg_region & # on reg region chrom
                                                                                    same_TAD_variant_and_reg_region & # in reg region TAD
                                                                                    (variant_position < reg_region_start), # low of reg region
                                                                                    "ON_REG_CHROM__IN_REG_TAD__LOW_OF_REG",
                                                                                    ifelse(
                                                                                        (!same_chrom_gene_and_reg_region) & (!same_chrom_gene_and_variant) & same_chrom_variant_and_reg_region & # on reg region chrom
                                                                                        (!same_TAD_variant_and_reg_region) & # outside reg region TAD
                                                                                        (variant_position > reg_region_end), # high of reg region
                                                                                        "ON_REG_CHROM__OUTSIDE_TAD__HIGH_OF_REG",
                                                                                        ifelse(
                                                                                            (!same_chrom_gene_and_reg_region) & (!same_chrom_gene_and_variant) & same_chrom_variant_and_reg_region & # on reg region chrom
                                                                                            (!same_TAD_variant_and_reg_region) & # outside reg region TAD
                                                                                            (variant_position < reg_region_start), # low of reg region
                                                                                            "ON_REG_CHROM__OUTSIDE_TAD__LOW_OF_REG",
                                                                                            ifelse(
                                                                                                (!same_chrom_gene_and_reg_region) & same_chrom_gene_and_variant & (!same_chrom_variant_and_reg_region) & # on gene chrom
                                                                                                same_TAD_gene_and_variant & # in gene TAD
                                                                                                (variant_position >= gene_start) & (variant_position <= gene_end),  # within gene
                                                                                                "ON_GENE_CHROM__IN_GENE_TAD__IN_GENE",
                                                                                                ifelse(
                                                                                                    (!same_chrom_gene_and_reg_region) & same_chrom_gene_and_variant & (!same_chrom_variant_and_reg_region) & # on gene chrom
                                                                                                    same_TAD_gene_and_variant & # in gene TAD
                                                                                                    (variant_position > gene_end),  # high of gene
                                                                                                    "ON_GENE_CHROM__IN_GENE_TAD__HIGH_OF_GENE",
                                                                                                    ifelse(
                                                                                                        (!same_chrom_gene_and_reg_region) & same_chrom_gene_and_variant & (!same_chrom_variant_and_reg_region) & # on gene chrom
                                                                                                        same_TAD_gene_and_variant & # in gene TAD
                                                                                                        (variant_position < gene_start),  # low of gene
                                                                                                        "ON_GENE_CHROM__IN_GENE_TAD__LOW_OF_GENE",
                                                                                                        ifelse(
                                                                                                            (!same_chrom_gene_and_reg_region) & same_chrom_gene_and_variant & (!same_chrom_variant_and_reg_region) & # on gene chrom
                                                                                                            (!same_TAD_gene_and_variant) & # outside gene TAD
                                                                                                            (variant_position > gene_end),  # high of gene
                                                                                                            "ON_GENE_CHROM__OUTSIDE_TAD__HIGH_OF_GENE",
                                                                                                            ifelse(
                                                                                                                (!same_chrom_gene_and_reg_region) & same_chrom_gene_and_variant & (!same_chrom_variant_and_reg_region) & # on gene chrom
                                                                                                                (!same_TAD_gene_and_variant) & # outside gene TAD
                                                                                                                (variant_position < gene_start),  # low of gene
                                                                                                                "ON_GENE_CHROM__OUTSIDE_TAD__LOW_OF_GENE",
                                                                                                                ifelse(
                                                                                                                    (!same_chrom_gene_and_variant) & (!same_chrom_variant_and_reg_region), # other chrom
                                                                                                                    "OTHER_CHROM",
                                                                                                                    "NOT_POSSIBLE"
                                                                                                                )
                                                                                                            )
                                                                                                        )
                                                                                                    )
                                                                                                )
                                                                                            )
                                                                                        )
                                                                                    )
                                                                                )
                                                                            )
                                                                        )
                                                                    )
                                                                )
                                                            )
                                                        )
                                                    )
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
    )

    return(position_category)
}

evaluate_CNA_genehancer_variant_position <- function(CNA_genehancer_combinations) {
    if(!(nrow(CNA_genehancer_combinations)==0)){
        GenomeInfoDb::seqlevelsStyle(CNA_genehancer_combinations$cna_chrom) <- "UCSC"
    }

    position_category_start <- evaluate_variant_position_by_single_location(
        variant_chrom = CNA_genehancer_combinations$cna_chrom,
        variant_position = CNA_genehancer_combinations$cna_start,
        gene_chrom = CNA_genehancer_combinations$chrom_gene,
        gene_start = CNA_genehancer_combinations$start_gene,
        gene_end = CNA_genehancer_combinations$end_gene,
        reg_region_chrom = CNA_genehancer_combinations$chrom_genehancer,
        reg_region_start = CNA_genehancer_combinations$start_genehancer,
        reg_region_end = CNA_genehancer_combinations$end_genehancer,
        gene_TAD_start = CNA_genehancer_combinations$min_TAD_start,
        gene_TAD_end = CNA_genehancer_combinations$max_TAD_end,
        reg_region_TAD_start = CNA_genehancer_combinations$min_TAD_start_genehancer,
        reg_region_TAD_end = CNA_genehancer_combinations$max_TAD_end_genehancer
    )

    position_category_end <- evaluate_variant_position_by_single_location(
        variant_chrom = CNA_genehancer_combinations$cna_chrom,
        variant_position = CNA_genehancer_combinations$cna_end,
        gene_chrom = CNA_genehancer_combinations$chrom_gene,
        gene_start = CNA_genehancer_combinations$start_gene,
        gene_end = CNA_genehancer_combinations$end_gene,
        reg_region_chrom = CNA_genehancer_combinations$chrom_genehancer,
        reg_region_start = CNA_genehancer_combinations$start_genehancer,
        reg_region_end = CNA_genehancer_combinations$end_genehancer,
        gene_TAD_start = CNA_genehancer_combinations$min_TAD_start,
        gene_TAD_end = CNA_genehancer_combinations$max_TAD_end,
        reg_region_TAD_start = CNA_genehancer_combinations$min_TAD_start_genehancer,
        reg_region_TAD_end = CNA_genehancer_combinations$max_TAD_end_genehancer
    )
    
    CNA_genehancer_combinations[, position_category_CNA := paste0(position_category_start,"->",position_category_end)]
}


evaluate_CNA_chipseq_variant_position <- function(CNA_chipseq_combinations) {
    if(!(nrow(CNA_chipseq_combinations)==0)){
        GenomeInfoDb::seqlevelsStyle(CNA_chipseq_combinations$cna_chrom) <- "UCSC"
    }

    position_category_start <- evaluate_variant_position_by_single_location(
        variant_chrom = CNA_chipseq_combinations$cna_chrom,
        variant_position = CNA_chipseq_combinations$cna_start,
        gene_chrom = CNA_chipseq_combinations$chrom,
        gene_start = CNA_chipseq_combinations$start,
        gene_end = CNA_chipseq_combinations$end,
        reg_region_chrom = CNA_chipseq_combinations$chipseq_chrom,
        reg_region_start = CNA_chipseq_combinations$chipseq_start,
        reg_region_end = CNA_chipseq_combinations$chipseq_end,
        gene_TAD_start = CNA_chipseq_combinations$min_TAD_start,
        gene_TAD_end = CNA_chipseq_combinations$max_TAD_end,
        reg_region_TAD_start = CNA_chipseq_combinations$min_TAD_start_chipseq,
        reg_region_TAD_end = CNA_chipseq_combinations$max_TAD_end_chipseq
    )

    position_category_end <- evaluate_variant_position_by_single_location(
        variant_chrom = CNA_chipseq_combinations$cna_chrom,
        variant_position = CNA_chipseq_combinations$cna_end,
        gene_chrom = CNA_chipseq_combinations$chrom,
        gene_start = CNA_chipseq_combinations$start,
        gene_end = CNA_chipseq_combinations$end,
        reg_region_chrom = CNA_chipseq_combinations$chipseq_chrom,
        reg_region_start = CNA_chipseq_combinations$chipseq_start,
        reg_region_end = CNA_chipseq_combinations$chipseq_end,
        gene_TAD_start = CNA_chipseq_combinations$min_TAD_start,
        gene_TAD_end = CNA_chipseq_combinations$max_TAD_end,
        reg_region_TAD_start = CNA_chipseq_combinations$min_TAD_start_chipseq,
        reg_region_TAD_end = CNA_chipseq_combinations$max_TAD_end_chipseq
    )
    
    CNA_chipseq_combinations[, position_category_CNA := paste0(position_category_start,"->",position_category_end)]
}


evaluate_SV_genehancer_variant_position <- function(SV_genehancer_combinations) {

    if(!(nrow(SV_genehancer_combinations)==0)){
        GenomeInfoDb::seqlevelsStyle(SV_genehancer_combinations$sv_break1_chrom) <- "UCSC"
        GenomeInfoDb::seqlevelsStyle(SV_genehancer_combinations$sv_break2_chrom) <- "UCSC"
    }

    position_category_break1 <- evaluate_variant_position_by_single_location(
        variant_chrom = SV_genehancer_combinations$sv_break1_chrom,
        variant_position = SV_genehancer_combinations$sv_break1_pos,
        gene_chrom = SV_genehancer_combinations$chrom_gene,
        gene_start = SV_genehancer_combinations$start_gene,
        gene_end = SV_genehancer_combinations$end_gene,
        reg_region_chrom = SV_genehancer_combinations$chrom_genehancer,
        reg_region_start = SV_genehancer_combinations$start_genehancer,
        reg_region_end = SV_genehancer_combinations$end_genehancer,
        gene_TAD_start = SV_genehancer_combinations$min_TAD_start,
        gene_TAD_end = SV_genehancer_combinations$max_TAD_end,
        reg_region_TAD_start = SV_genehancer_combinations$min_TAD_start_genehancer,
        reg_region_TAD_end = SV_genehancer_combinations$max_TAD_end_genehancer
    )

    position_category_break2 <- evaluate_variant_position_by_single_location(
        variant_chrom = SV_genehancer_combinations$sv_break2_chrom,
        variant_position = SV_genehancer_combinations$sv_break2_pos,
        gene_chrom = SV_genehancer_combinations$chrom_gene,
        gene_start = SV_genehancer_combinations$start_gene,
        gene_end = SV_genehancer_combinations$end_gene,
        reg_region_chrom = SV_genehancer_combinations$chrom_genehancer,
        reg_region_start = SV_genehancer_combinations$start_genehancer,
        reg_region_end = SV_genehancer_combinations$end_genehancer,
        gene_TAD_start = SV_genehancer_combinations$min_TAD_start,
        gene_TAD_end = SV_genehancer_combinations$max_TAD_end,
        reg_region_TAD_start = SV_genehancer_combinations$min_TAD_start_genehancer,
        reg_region_TAD_end = SV_genehancer_combinations$max_TAD_end_genehancer
    )

    SV_genehancer_combinations[, position_category_SV := paste0(position_category_break1,"<->",position_category_break2)]
    return(SV_genehancer_combinations)
}


evaluate_SV_chipseq_variant_position <- function(SV_chipseq_combinations) {

    if(!(nrow(SV_chipseq_combinations)==0)){
        GenomeInfoDb::seqlevelsStyle(SV_chipseq_combinations$sv_break1_chrom) <- "UCSC"
        GenomeInfoDb::seqlevelsStyle(SV_chipseq_combinations$sv_break2_chrom) <- "UCSC"
    }

    position_category_break1 <- evaluate_variant_position_by_single_location(
        variant_chrom = SV_chipseq_combinations$sv_break1_chrom,
        variant_position = SV_chipseq_combinations$sv_break1_pos,
        gene_chrom = SV_chipseq_combinations$chrom,
        gene_start = SV_chipseq_combinations$start,
        gene_end = SV_chipseq_combinations$end,
        reg_region_chrom = SV_chipseq_combinations$chipseq_chrom,
        reg_region_start = SV_chipseq_combinations$chipseq_start,
        reg_region_end = SV_chipseq_combinations$chipseq_end,
        gene_TAD_start = SV_chipseq_combinations$min_TAD_start,
        gene_TAD_end = SV_chipseq_combinations$max_TAD_end,
        reg_region_TAD_start = SV_chipseq_combinations$min_TAD_start_chipseq,
        reg_region_TAD_end = SV_chipseq_combinations$max_TAD_end_chipseq
    )

    position_category_break2 <- evaluate_variant_position_by_single_location(
        variant_chrom = SV_chipseq_combinations$sv_break2_chrom,
        variant_position = SV_chipseq_combinations$sv_break2_pos,
        gene_chrom = SV_chipseq_combinations$chrom,
        gene_start = SV_chipseq_combinations$start,
        gene_end = SV_chipseq_combinations$end,
        reg_region_chrom = SV_chipseq_combinations$chipseq_chrom,
        reg_region_start = SV_chipseq_combinations$chipseq_start,
        reg_region_end = SV_chipseq_combinations$chipseq_end,
        gene_TAD_start = SV_chipseq_combinations$min_TAD_start,
        gene_TAD_end = SV_chipseq_combinations$max_TAD_end,
        reg_region_TAD_start = SV_chipseq_combinations$min_TAD_start_chipseq,
        reg_region_TAD_end = SV_chipseq_combinations$max_TAD_end_chipseq
    )


    

    SV_chipseq_combinations[, position_category_SV := paste0(position_category_break1,"<->",position_category_break2)]
    return(SV_chipseq_combinations)
}