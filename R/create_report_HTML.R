create_report_HTML <- function(processed_data,
                               toggle_plot_HTML,
                               SV_recurrent_TAD_combinations_HTML,
                               HTML_mut_analysis,
                               analysis_by_gene_HTML) {
    paste0(
        '
        <html>
            <head>
                <title>Revana Report</title>
                <link
                    href="https://cdn.jsdelivr.net/npm/bootstrap@5.0.1/dist/css/bootstrap.min.css"
                    rel="stylesheet"
                    integrity="sha384-+0n0xVW2eSR5OomGNYDnhzAbDsOXxcvSN1TPprVMTNDbiYZCxYbOOl7+AMvyTG2x"
                    crossorigin="anonymous"
                />
                <script
                    src="https://cdn.jsdelivr.net/npm/bootstrap@5.0.1/dist/js/bootstrap.bundle.min.js"
                    integrity="sha384-gtEjrD/SeCtmISkJkNUaaKMoLD0//ElJ19smozuHV6z3Iehds+3Ulb9Bn9Plx0x4"
                    crossorigin="anonymous"
                ></script>
                <style>
                    .sv-pos-table,
                    .cna-pos-table {
                        table-layout: fixed;
                    }

                    .sv-pos-table td:first-child,
                    .sv-pos-table th:first-child  {
                        width: 15em;
                    }

                    picture { display: none; }

                    @media (min-width: 50px) {
                        picture {
                        display: block;
                        }
                    }
                </style>
            </head>
            <body>
                <div class="container" style="max-width: 55em;">
                    <div class="vh-100 d-flex justify-content-center align-items-center">
                        <div class="bg-white p-5 rounded shadow">
                            <h1 class="display-3">Revana Analysis Report</h1>
                            <h5 class="display-5">', format(Sys.Date(), "%Y-%m-%d"), '</h5>
                            <p>Revana was developed by<br>Elias Ulrich, Stefan Pfister and Natalie Jaeger</p>
                        </div>
                    </div>

                    <h1 class="mb-5">Part I: Input data</h1>

                    <div class = "card my-4">
                        <h5 class="card-header">Analysed samples
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 text-center">
                            <div class="accordion accordion-flush" id="accordion-analysed-samples">
                                <div class="accordion-item">
                                    <div class="accordion-header" id="accordion-heading-number-of-samples">
                                        <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapse-sample-ids">
                                            <div class="flex-grow-1 text-center">
                                                <h5>Number of analysed samples</h5>
                                                <h3>',
        length(unique(processed_data$result_paths$sample_id)),
        '</h3>
                                            </div>
                                        </button>
                                    </div>
                                    <div id="collapse-sample-ids" class="accordion-collapse collapse"  style = "max-height: 20rem; overflow-y: scroll;" data-bs-parent="#accordion-analysed-samples">
                                        <ul class="list-group list-group-flush">',
        paste0("<li class=\"list-group-item\">", processed_data$result_paths$sample_id, "</li>", collapse = "\n"),
        '</ul>
                                    </div>
                                </div>
                                <div class="accordion-item">
                                    <div class="accordion-header" id="accordion-heading-number-of-processed_data$subgroups">
                                        <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapse-subgroup-names">
                                            <div class="flex-grow-1 text-center">
                                                <h5>Number of analysed subgroups</h5>
                                                <h3>',
        length(processed_data$subgroups_table$subgroup),
        '</h3>
                                            </div>
                                        </button>
                                    </div>
                                    <div id="collapse-subgroup-names" class="accordion-collapse collapse"  style = "max-height: 20rem; overflow-y: scroll;" data-bs-parent="#accordion-analysed-samples">
                                        <ul class="list-group list-group-flush">',
        paste0("<li class=\"list-group-item\">", processed_data$subgroups_table$subgroup, " - ", processed_data$subgroups_table$n, " samples</li>", collapse = "\n"),
        '</ul>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>

                    <div class = "card my-4">
                        <h5 class="card-header">Sample Subgroup Composition
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 text-center">
                            <img class ="w-100" src="figures/subgroups_bar_chart.png" />
                        </div>
                    </div>

                    <div class = "card my-4">
                        <h5 class="card-header">SNP Markers
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 border-bottom">
                          SNP marker input contains all somatic and germline single nucleotide polymorphisms of the tumor sample. It is used to calculate allele-specific expression.
                        </div>
                        <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<img class =\"w-50 p-2 \" src=\"figures/marker_bar_chart_bygroup_", processed_data$subgroups, ".png\" />", collapse = "\n"),
        '</div>
                    </div>

                    <div class = "card my-4">
                        <h5 class="card-header">Somatic SNVs/InDels
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 border-bottom">
                          Somatic SNV input contains somatic single nucleotide variants (SNVs) and small insertion and deletions (InDels) of the tumor sample. These variants are being evaluated for being candidate regulatory variants.
                        </div>
                        <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<img class =\"w-50 p-2 \" src=\"figures/somatic_SNV_bar_chart_bygroup_", processed_data$subgroups, ".png\" />", collapse = "\n"),
        '</div>
                    </div>

                    <div class = "card my-4">
                        <h5 class="card-header">SVs
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 border-bottom">
                          SV input contains somatic structural variants of the tumor sample. These variants are being evaluated for being candidate regulatory variants. They are also tested for recurrent SV TAD juxtaposition.
                        </div>
                        <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<img class =\"w-50 p-2 \" src=\"figures/SV_bar_chart_bygroup_", processed_data$subgroups, ".png\" />", collapse = "\n"),
        '
                        </div>
                    </div>

                    <div class = "card my-4">
                        <h5 class="card-header">CNAs
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 border-bottom">
                          CNA input contains somatic copy number alterations of the tumor sample. These variants are being evaluated for being candidate regulatory variants.
                        </div>
                        <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<img class =\"w-50 p-2 \" src=\"figures/CNA_bar_chart_bygroup_", processed_data$subgroups, ".png\" />", collapse = "\n"),
        '
                        </div>
                    </div>

                    <div class = "card my-4">
                        <h5 class="card-header">Gene Expression
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 border-bottom">
                          Gene Expression input contains RNA-Seq feature counts in FPKM. Use of this input includes Outlier High Expression (OHE) detection.
                        </div>
                        <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<img class =\"w-50 p-2 \" src=\"figures/expression_box_plot_bygroup_", processed_data$subgroups, ".png\" />", collapse = "\n"),
        '
                        </div>
                    </div>

                    <div class = "card my-4">
                        <h5 class="card-header">Gene Expression - PCA
                            <span class="badge bg-primary float-end">Input data</span>
                        </h5>
                        <div class="p-3 border-bottom">
                          These graphs show the principal component analysis (PCA) of gene expression of the included subgroups. They help to quickly detect batch effect distortions between subsets of the samples, and thus to assert comparability between the samples.
                        </div>
                        <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<img class =\"w-50 p-2 \" src=\"figures/expression_PCA_plot_bygroup_", processed_data$subgroups, ".png\" />", collapse = "\n"),
        '
                        </div>
                    </div>

                    <h1 class="mb-5">Part II: Results Overview</h1>

                    <div class = "card my-4">
                        <h5 class="card-header">Cis-Activation
                            <span class="badge bg-cyan float-end">Overview</span>
                        </h5>
                        <div class="p-3 text-center">
                            <div class="accordion accordion-flush" id="accordion-samples-with-cis-activation">
                                <div class="accordion-item">
                                    <div class="accordion-header" id="accordion-heading-number-of-samples-with-cis-activation">
                                        <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapse-sample-ids-with-cis-activation">
                                            <div class="flex-grow-1 text-center">
                                                <h5>Samples <strong>with</strong> detected cis-activation</h5>
                                                <h3>', length((processed_data$cis_activated_genes_by_sample %>% dplyr::filter(n_cis_activated_genes > 0))$sample_ID), '</h3>
                                            </div>
                                        </button>
                                    </div>
                                    <div id="collapse-sample-ids-with-cis-activation" class="accordion-collapse collapse"  style = "max-height: 20rem; overflow-y: scroll;" data-bs-parent="#accordion-samples-with-cis-activation">
                                        <ul class="list-group list-group-flush">',
        paste0("<li class=\"list-group-item\">", (processed_data$cis_activated_genes_by_sample %>% dplyr::filter(n_cis_activated_genes > 0))$sample_ID, "</li>", collapse = "\n"),
        '</ul>
                                    </div>
                                </div>
                                <div class="accordion-item">
                                    <div class="accordion-header" id="accordion-heading-number-of-samples-without-cis-activation">
                                        <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" data-bs-target="#collapse-sample-ids-without-cis-activation">
                                            <div class="flex-grow-1 text-center">
                                                <h5>Samples <strong>without</strong> detected cis-activation</h5>
                                                <h3>', length((processed_data$cis_activated_genes_by_sample %>% dplyr::filter(n_cis_activated_genes == 0))$sample_ID), '</h3>
                                            </div>
                                        </button>
                                    </div>
                                    <div id="collapse-sample-ids-without-cis-activation" class="accordion-collapse collapse"  style = "max-height: 20rem; overflow-y: scroll;" data-bs-parent="#accordion-samples-with-cis-activation">
                                        <ul class="list-group list-group-flush">',
        paste0("<li class=\"list-group-item\">", (processed_data$cis_activated_genes_by_sample %>% dplyr::filter(n_cis_activated_genes == 0))$sample_ID, "</li>", collapse = "\n"),
        '</ul>
                                    </div>
                                </div>
                                <div class="text-center p-3 border-bottom">
                                    <h5>Cis-activated genes</h5>
                                    <h3>', sum(processed_data$cis_activation_summary_table$cis_activated_gene, na.rm = TRUE), '</h3>
                                </div>
                                <div class="text-center p-3 border-bottom">
                                    <h5>Unique cis-activated genes</h5>
                                    <h3>', length(unique(processed_data$cis_activation_summary_table$gene_name[(!is.na(processed_data$cis_activation_summary_table$cis_activated_gene)) & processed_data$cis_activation_summary_table$cis_activated_gene])), '</h3>
                                </div>
                                <div class="text-center p-3 border-bottom">
                                    <h5>Average cis-activated genes per sample</h5>
                                    <h3>', sum(processed_data$cis_activation_summary_table$cis_activated_gene, na.rm = TRUE) / length(unique(processed_data$result_paths$sample_id)), '</h3>
                                </div>
                            </div>
                        </div>
                        <div class="p-3 text-center  d-flex flex-wrap">',
        paste0("<img class =\"w-50 p-2 \" src=\"figures/cis_activated_genes_bar_chart_bygroup_", processed_data$subgroups, ".png\" />", collapse = "\n"),
        '
                        </div>
                    </div>

                    <div class = "card my-4">
                      <h5 class="card-header">Algorithm insights
                          <span class="badge bg-danger float-end">Results</span>
                      </h5>
                      <div class="p-3 text-center border-bottom">
                        <button
                          id="show-algorithm-insights-button"
                          class="btn btn-secondary btn-js-show btn-js-enable btn-js-disable-self"
                          data-enable-id="hide-algorithm-insights-button"
                          data-show-id="algorithm-insights"
                        >Show</button>
                        <button
                          id="hide-algorithm-insights-button"
                          class="btn btn-secondary btn-js-hide btn-js-enable btn-js-disable-self"
                          data-enable-id="show-algorithm-insights-button"
                          data-hide-id="algorithm-insights" disabled="true"
                        >Hide</button>
                      </div>
                      <div id="algorithm-insights" class="p-3 bg-secondary" style="display:none;">

                        <div class = "card my-4">
                            <h5 class="card-header">Type of ASE detection
                                <span class="badge bg-danger float-end">Results</span>
                            </h5>
                            <div class="p-3 text-center border-bottom">
                                <p>Type of ASE detection - Whole Cohort</p>
                                <img class ="w-100 p-2" src="figures/ase_detection_upset_plot_whole_cohort.svg" />
                            </div>
                            <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<div class =\"w-50 p-2 \"><p>Type of ASE detection - ", processed_data$subgroups, "</p><img class =\"w-100\" src=\"figures/ase_detection_upset_plot_bygroup_", processed_data$subgroups, ".svg\" /></div>", collapse = "\n"),
        '
                            </div>
                        </div>

                        <div class = "card my-4">
                            <h5 class="card-header">Applied Filters
                                <span class="badge bg-danger float-end">Results</span>
                            </h5>
                            <div class="p-3 text-center border-bottom">
                                <p>Applied Filters - Whole Cohort</p>
                                <img class ="w-100 p-2" src="figures/applied_filters_upset_plot_whole_cohort.svg" />
                            </div>
                            <div class="p-3 text-center d-flex flex-wrap">',
        paste0("<div class =\"w-50 p-2 \"><p>Applied Filters - ", processed_data$subgroups, "</p><img class =\"w-100\" src=\"figures/applied_filters_upset_plot_bygroup_", processed_data$subgroups, ".svg\" /></div>", collapse = "\n"),
        '
                            </div>
                        </div>



                        <div class = "card my-4">
                            <h5 class="card-header">ASE and OHE values for cis-activated genes
                                <span class="badge bg-danger float-end">Results</span>
                            </h5>
                            <div class="p-3 text-center">
                                <img class ="w-100" src="figures/OHE_ASE_q_value_dot_plot.png" />
                            </div>
                        </div>
                      </div>
                    </div>

                    ', toggle_plot_HTML, SV_recurrent_TAD_combinations_HTML, '


                    <h1 class="mb-5">Part III: Analysis of mutation types</h1>
                    ', HTML_mut_analysis, '

                    <h1 class="mb-5">Part IV: Analysis by gene</h1>
                    <div class="mb-4 p-3 border border-secondary rounded">
                        <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class="analysis_by_gene_not_proteine_coding"
                            data-hide-class=""
                            data-enable-id="btn-hide-analysis_by_gene_not_proteine_coding"
                            id="btn-show-analysis_by_gene_not_proteine_coding" disabled="true">All Gene Types</button>
                        <button
                            class="btn btn-secondary btn-js-show btn-js-hide btn-js-disable-self btn-js-enable"
                            data-show-class=""
                            data-hide-class="analysis_by_gene_not_proteine_coding"
                            data-enable-id="btn-show-analysis_by_gene_not_proteine_coding"
                            id="btn-hide-analysis_by_gene_not_proteine_coding">Protein coding genes only</button>
                    </div>
                    ',
        analysis_by_gene_HTML,
        '
                </div>
                <script>
                    function jsShow(e) {
                    const clickedElement = e.target;
                    const dataShowId = clickedElement.getAttribute("data-show-id");
                    const dataShowClass = clickedElement.getAttribute("data-show-class");
                    if(dataShowId){
                        const elementWithId = document.getElementById(dataShowId);
                        elementWithId.style.display = "block";
                    }
                    if(dataShowClass){
                        const elementsWithClass = document.getElementsByClassName(dataShowClass);
                        if (elementsWithClass.length) {
                        for (let elementWithClass of elementsWithClass) {
                            elementWithClass.style.display = "block";
                        }
                        }
                    }
                    }
                    function jsHide(e) {
                    const clickedElement = e.target;
                    const dataHideId = clickedElement.getAttribute("data-hide-id");
                    const dataHideClass = clickedElement.getAttribute("data-hide-class");
                    if(dataHideId){
                        const elementWithId = document.getElementById(dataHideId);
                        elementWithId.style.display = "none";
                    }
                    if(dataHideClass){
                        const elementsWithClass = document.getElementsByClassName(dataHideClass);
                        if (elementsWithClass.length) {
                        for (let elementWithClass of elementsWithClass) {
                            elementWithClass.style.display = "none";
                        }
                        }
                    }
                    }

                    function jsDisable(e) {
                    const clickedElement = e.target;
                    const dataDisableId = clickedElement.getAttribute("data-disable-id");
                    const dataDisableClass = clickedElement.getAttribute("data-disable-class");
                    if(dataDisableId){
                        const elementWithId = document.getElementById(dataDisableId);
                        elementWithId.setAttribute("disabled","true");
                    }
                    if(dataDisableClass){
                        const elementsWithClass = document.getElementsByClassName(dataDisableClass);
                        if (elementsWithClass.length) {
                        for (let elementWithClass of elementsWithClass) {
                            elementWithClass.setAttribute("disabled","true");
                        }
                        }
                    }
                    }

                    function jsEnable(e) {
                    const clickedElement = e.target;
                    const dataEnableId = clickedElement.getAttribute("data-enable-id");
                    const dataEnableClass = clickedElement.getAttribute("data-enable-class");
                    if(dataEnableId){
                        const elementWithId = document.getElementById(dataEnableId);
                        elementWithId.removeAttribute("disabled");
                    }
                    if(dataEnableClass){
                        const elementsWithClass = document.getElementsByClassName(dataEnableClass);
                        if (elementsWithClass.length) {
                        for (let elementWithClass of elementsWithClass) {
                            elementWithClass.removeAttribute("disabled");
                        }
                        }
                    }
                    }

                    function jsDisableSelf(e) {
                    const clickedElement = e.target;
                    clickedElement.setAttribute("disabled","true");
                    }

                    const jsShowBtns = document.getElementsByClassName("btn-js-show");
                    const jsHideBtns = document.getElementsByClassName("btn-js-hide");
                    const jsDisableSelfBtns = document.getElementsByClassName("btn-js-disable-self");
                    const jsDisableBtns = document.getElementsByClassName("btn-js-disable");
                    const jsEnableBtns = document.getElementsByClassName("btn-js-enable");

                    for (let jsShowBtn of jsShowBtns) {
                    jsShowBtn.addEventListener("click",jsShow)
                    }
                    for (let jsHideBtn of jsHideBtns) {
                    jsHideBtn.addEventListener("click",jsHide)
                    }
                    for (let jsDisableSelfBtn of jsDisableSelfBtns) {
                    jsDisableSelfBtn.addEventListener("click",jsDisableSelf)
                    }
                    for (let jsDisableBtn of jsDisableBtns) {
                    jsDisableBtn.addEventListener("click",jsDisable)
                    }
                    for (let jsEnableBtn of jsEnableBtns) {
                    jsEnableBtn.addEventListener("click",jsEnable)
                    }
                </script>
            </body>
        </html>
    '
    )
}