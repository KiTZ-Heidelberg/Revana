create_HTML_for_SV_recurrent_TAD_combinations <- function(processed_data) {
    paste0(
        '<div class = "card my-4">
            <h5 class="card-header">SVs: Recurrent SV TAD juxtaposition
                <span class="badge bg-danger float-end">Results</span>
            </h5>
            <div class="card-body p-3 border-bottom">
                <div style="display: block;
                        overflow-y: scroll;
                        height: 30em;">
                    ', knitr::kable(
            dplyr::select(processed_data$recurrent_SV_TAD_combinations_summary, TADs = TAD_combination, `N of samples` = n_samples_with_TAD_combination, Samples = samples_with_TAD_combination, `N of samples with cis-activated genes affected` = n_samples_with_TAD_combination_and_cis_activated_genes_affected, `cis-activated genes` = affected_cis_activated_genes),
            format = "html",
            table.attr = "class=\"table table-striped table-bordered\""
        ), "
                </div>
        </div>
        </div>"
    )
}