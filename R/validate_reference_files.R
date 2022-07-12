validate_reference_files <- function(gene_annotation_ref_file_path,
                                     gene_annotation_exons_ref_file_path,
                                     fimo_motif_ref_path,
                                     motif_id_tf_gene_name_table_path,
                                     TAD_file_path,
                                     genehancer_ref_file_path,
                                     chipseq_file_path,
                                     run_tf_binding_site_analysis) {
    # check for existence ----------------------
    # gene annotation reference file
    check_file_existence(gene_annotation_ref_file_path, name_of_file_type = "gene annotation reference file")

    # gene annotation exons reference file
    check_file_existence(gene_annotation_exons_ref_file_path, name_of_file_type = "gene annotation exons reference file")

    if (run_tf_binding_site_analysis == TRUE) {
        # fimo motif reference file
        check_file_existence(fimo_motif_ref_path, name_of_file_type = "fimo motif reference file")

        # motif id tf gene name table
        check_file_existence(motif_id_tf_gene_name_table_path, name_of_file_type = "motif_id tf_gene_name table reference file")
    }

    # TAD file
    check_file_existence(TAD_file_path, name_of_file_type = "TAD file")

    # genehancer reference file
    if (is.null(genehancer_ref_file_path)) {
        cat("Revana is running WITHOUT GeneHancer functionality...\n")
    } else {
        check_file_existence(genehancer_ref_file_path, name_of_file_type = "genehancer reference file")
    }

    # chipseq file
    if (is.null(chipseq_file_path)) {
        cat("Revana is running WITHOUT ChIP-Seq functionality...\n")
    } else {
        check_file_existence(chipseq_file_path, name_of_file_type = "ChIP-Seq file")
    }


    # check for right headers ------------------
    # gene annotation reference file
    required_cols_gene_annotation_ref_file <- c("chrom", "start", "end", "width", "strand", "gene_name", "gene_type", "imprinting_status", "imprinting_expressed_allele", "cancer_gene_role_in_cancer", "is_cancer_gene")
    check_file_table_header(gene_annotation_ref_file_path, required_cols = required_cols_gene_annotation_ref_file, name_of_file_type = "gene annotation reference file")

    # gene annotation reference exons file
    required_cols_gene_annotation_exons_ref_file <- c("chrom", "start", "end", "gene_name")
    check_file_table_header(gene_annotation_exons_ref_file_path, required_cols = required_cols_gene_annotation_exons_ref_file, name_of_file_type = "gene annotation exons reference file")

    # motif id tf gene name table
    if (run_tf_binding_site_analysis == TRUE) {
        required_cols_motif_id_tf_gene_name_table <- c("motif_id", "tf_gene_name")
        check_file_table_header(motif_id_tf_gene_name_table_path, required_cols = required_cols_motif_id_tf_gene_name_table, name_of_file_type = "motif_id tf_gene_name table reference file")
    }


    # TAD file
    required_cols_TAD_file <- c("chrom", "start", "end")
    check_file_table_header(TAD_file_path, required_cols = required_cols_TAD_file, name_of_file_type = "TAD file")

    # genehancer reference file
    if (!is.null(genehancer_ref_file_path)) {
        required_cols_genehancer_ref_file <- c("chrom", "feature_name", "start", "end", "score", "genehancer_id", "connected_gene", "connected_gene_score", "is_elite", "is_association_elite")
        check_file_table_header(genehancer_ref_file_path, required_cols = required_cols_genehancer_ref_file, name_of_file_type = "genehancer reference file")
    }


    # chipseq reference file
    if (!is.null(chipseq_file_path)) {
        required_cols_chipseq_file <- c("chrom", "start", "end", "cluster")
        check_file_table_header(chipseq_file_path, required_cols = required_cols_chipseq_file, name_of_file_type = "ChIP-Seq file")
    }
}

validate_gene_annotation_ref_file_path <- function(gene_annotation_ref_file_path) {
    check_file_existence(gene_annotation_ref_file_path, name_of_file_type = "Gene Annotation Reference File")
    required_cols_gene_annotation_ref_file <- c("chrom", "start", "end", "width", "strand", "gene_name", "gene_type", "imprinting_status", "imprinting_expressed_allele", "cancer_gene_role_in_cancer", "is_cancer_gene")
    check_file_table_header(gene_annotation_ref_file_path, required_cols = required_cols_gene_annotation_ref_file, name_of_file_type = "Gene Annotation Reference File")
}