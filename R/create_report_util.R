format_number <- function(x) {
    return(format(x, big.mark = " "))
}


write_HTML_to_file <- function(html_content,
                               write_path,
                               print_before = "Writing HTML to file...",
                               print_after = "NEW REPORT CREATED \n") {
    cat(print_before)

    con <- file(write_path, "w")
    writeLines(html_content, con)
    close(con)

    cat(print_after)
}

create_dir_if_missing <- function(path) {
    if (!dir.exists(path)) {
        dir.create(path)
    }
    return(path)
}

create_file_if_missing <- function(path) {
    if (!file.exists(path)) {
        file.create(path)
    }
    return(path)
}

get_n_colors <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

create_color_palette <- function(subgroups) {
    subgroup_colors <- get_n_colors(length(unique(subgroups)))
    names(subgroup_colors) <- unique(subgroups)
    return(
        list(
            scale_fill_subgroup = ggplot2::scale_fill_manual(name = "Subgroup", values = subgroup_colors),
            scale_color_subgroup = ggplot2::scale_color_manual(name = "Subgroup", values = subgroup_colors),
            scale_color_mutation_type = ggplot2::scale_color_manual(values = c(SV = "orange", somatic_SNV = "cyan"))
        )
    )
}

save_ggplot <- function(plot, filename) {
    ggplot2::ggsave(filename = filename, plot = plot)
}

cont_table <- function(a, b) {
    return(table(
        factor(a, levels = c(FALSE, TRUE)),
        factor(b, levels = c(FALSE, TRUE))
    ))
}