#' Dot plot showing top differential signaling LR-pairs
#' @description
#' The size of dots reflects significance of the LR-pairs, and the color
#' indicates the direction of the interaction.
#' @param dealr a `dealr` object. DO NOT FILTER.
#' @param sender_use Character vector of sender groups to be considered. Default
#' `NULL` includes all possible senders.
#' @param receiver_use Character vector of receiver groups to be considered.
#' Default `NULL` includes all possible receivers.
#' @param focus A single character string indicating a group of interest. If
#' provided, only show the communication involving this group and the plot
#' will be split to two panels for when this group is sending the ligand and
#' when this group is receiving the ligand. Default `NULL` use all possible
#' pairs of sender and receiver as provided with `sender_use` and
#' `receiver_use`.
#' @param pathway_use Character vector of pathways to filter the LR-pairs.
#' Applied after p-value filtering but prior to top n slicing. Default
#' `NULL`.
#' @param p_thresh a numeric value indicating the p-value threshold for
#' filtering the LR-pairs. Default is `0.01`.
#' @param top_n a numeric value indicating the number of top LR-pairs to be
#' plotted, ranked by p-value. Default `NULL` means all LR-pairs, usually
#' too many.
#' @inheritParams .dotplot
#' @inheritParams .theme_text_setter
#' @return ggplot object of the dot plot.
#' @export
plotLRPairDot <- function(
        dealr,
        sender_use = NULL,
        receiver_use = NULL,
        focus = NULL,
        pathway_use = NULL,
        p_thresh = 0.01,
        top_n = NULL,
        downreg_col = '#2166AC',
        upreg_col = '#B2182B',
        dot_size_range = c(0.3, 4),
        text_x_size = 8,
        text_y_size = 6,
        text_strip_size = 6,
        text_title_size = 10,
        text_legend_size = 8,
        text_legend_title_size = 8
) {
    if (!inherits(dealr, 'dealr')) {
        cli_abort('{.field dealr} must be a {.cls dealr} object')
    }
    sender_use <- sender_use %||% levels(dealr$sender)
    sender_use <- arg_match(sender_use, levels(dealr$sender), multiple = TRUE)
    receiver_use <- receiver_use %||% levels(dealr$receiver)
    receiver_use <- arg_match(receiver_use, levels(dealr$receiver), multiple = TRUE)

    dealr <- dealr %>% filter(
        .data[['sender']] %in% sender_use,
        .data[['receiver']] %in% receiver_use,
        .data[['p']] < p_thresh
    ) %>%
        mutate(
            pair = paste(.data[['sender']], '->', .data[['receiver']])
        )
    if (!is.null(focus)) {
        focus <- arg_match(focus, c(sender_use, receiver_use))
        dealr <- dealr %>%
            filter(.data[['sender']] == focus | .data[['receiver']] == focus) %>%
            mutate(
                role = ifelse(
                    .data[['sender']] == focus,
                    paste(focus, 'sending'),
                    paste(focus, 'receiving')
                )
            )
    }

    pathway_use <- pathway_use %||% unique(dealr[['pathway']])
    pathway_use <- arg_match(
        pathway_use, values = unique(dealr[['pathway']]), multiple = TRUE
    )
    dealr <- dealr %>% filter(.data[['pathway']] %in% pathway_use)

    if (!is.null(top_n)) {
        dealr <- dealr %>%
            ungroup() %>%
            arrange(.data[['p']]) %>%
            slice_head(n = top_n)
    }
    if (nrow(dealr) == 0) {
        cli_abort('No significant LR-pairs found for the sender {sender_use}.')
    }
    dealr %>%
        .dotplot(
            x = dealr$pair,
            y = dealr$LR_pair,
            size = -log10(dealr$p),
            fill = dealr$stat,
            downreg_col = downreg_col,
            upreg_col = upreg_col,
            fill_title = 'Z-score',
            size_title = '-log10(p-value)',
            dot_size_range = dot_size_range,
            x_position = 'bottom'
        ) +
        facet_grid(
            cols = (if ('role' %in% colnames(dealr)) vars(dealr$role) else NULL),
            rows = vars(dealr$pathway),
            scales = 'free', space = 'free_y'
        ) +
        theme(
            panel.spacing.y = unit(0, "lines"),
            strip.text.y.right = element_text(angle = 0, color = 'black'),
            strip.text.x.top = element_text(angle = 0, color = 'black'),
            strip.background = element_rect(fill = 'grey90', color = 'black', linewidth = 0.5),
        ) +
        .theme_text_setter(
            text_x_size = text_x_size,
            text_y_size = text_y_size,
            text_strip_size = text_strip_size,
            text_title_size = text_title_size,
            text_legend_size = text_legend_size,
            text_legend_title_size = text_legend_title_size
        )
}

#' Dot plot showing pathway enrichment from differential signaling LR-pairs
#' @description
#' The size of dots shows the fold enrichment of the pathways, and the color
#' indicates the significance and direction of the regulation of the pathways.
#' @param dealr_pe a `dealr_pe` object. DO NOT FILTER.
#' @param sender_use Character vector of sender groups to be considered. Default
#' `NULL` includes all possible senders.
#' @param receiver_use Character vector of receiver groups to be considered.
#' Default `NULL` includes all possible receivers.
#' @param focus A single character string indicating a group of interest. If
#' provided, only show the communication involving this group and the plot
#' will be split to two panels for when this group is sending the ligand and
#' when this group is receiving the ligand. Default `NULL` use all possible
#' pairs of sender and receiver as provided with `sender_use` and
#' `receiver_use`.
#' @param p_thresh a numeric value indicating the p-value threshold for
#' filtering the pathways. Default `0.01`.
#' @param top_n a numeric value indicating the number of top pathways to be
#' plotted, ranked by p-value. Default `NULL` means all pathways, usually
#' not too many.
#' @inheritParams .dotplot
#' @inheritParams .theme_text_setter
#' @return ggplot object of the dot plot.
#' @rdname plotDealrPE
#' @export
plotPathwayEnrichDot <- function(
        dealr_pe,
        sender_use = NULL,
        receiver_use = NULL,
        focus = NULL,
        p_thresh = 0.01,
        top_n = NULL,
        downreg_col = '#2166AC',
        upreg_col = '#B2182B',
        dot_size_range = c(0.3, 4),
        text_x_size = 8,
        text_y_size = 6,
        text_title_size = 10,
        text_legend_size = 8,
        text_legend_title_size = 8,
        text_strip_size = 6
) {
    if (!inherits(dealr_pe, 'dealr_pe')) {
        cli_abort('{.field dealr_pe} must be a {.cls dealr_pe} object')
    }
    sender_use <- sender_use %||% levels(dealr_pe$sender)
    sender_use <- arg_match(sender_use, levels(dealr_pe$sender), multiple = TRUE)
    receiver_use <- receiver_use %||% levels(dealr_pe$receiver)
    receiver_use <- arg_match(receiver_use, levels(dealr_pe$receiver), multiple = TRUE)

    dealr_pe <- dealr_pe %>% filter(
        .data[['sender']] %in% sender_use,
        .data[['receiver']] %in% receiver_use,
        .data[['p']] < p_thresh
    ) %>%
        mutate(
            pair = paste(.data[['sender']], '->', .data[['receiver']])
        )
    if (!is.null(focus)) {
        focus <- arg_match(focus, c(sender_use, receiver_use))
        dealr_pe <- dealr_pe %>%
            filter(.data[['sender']] == focus | .data[['receiver']] == focus) %>%
            mutate(
                role = ifelse(
                    .data[['sender']] == focus,
                    paste(focus, 'sending'),
                    paste(focus, 'receiving')
                )
            )
    }
    if (!is.null(top_n)) {
        dealr_pe <- dealr_pe %>%
            ungroup() %>%
            arrange(.data[['p']]) %>%
            slice_head(n = top_n)
    }
    if (nrow(dealr_pe) == 0) {
        cli_abort('No significant LR-pairs found for the sender {sender_use}.')
    }

    dealr_pe %>%
        .dotplot(
            x = dealr_pe$pair,
            y = dealr_pe$pathway,
            size = log(dealr_pe$fold_enrichment),
            fill = -log10(dealr_pe$p)*sign(dealr_pe$stat),
            downreg_col = downreg_col,
            upreg_col = upreg_col,
            fill_title = '-log10(p) x direction',
            size_title = 'log(Fold Enrichment)',
            dot_size_range = dot_size_range,
            x_position = 'bottom'
        ) +
        (
            if (!is.null(focus)) facet_grid(cols = vars(dealr_pe$role), space = 'free_x', scales = 'free_x')  else NULL
        ) +
        theme(
            panel.spacing.x = unit(0, "lines"),
            strip.text.x.top = element_text(angle = 0, color = 'black'),
            strip.background = element_rect(fill = 'grey90', color = 'black', linewidth = 0.5),
        ) +
        .theme_text_setter(
            text_x_size = text_x_size,
            text_y_size = text_y_size,
            text_title_size = text_title_size,
            text_legend_size = text_legend_size,
            text_legend_title_size = text_legend_title_size,
            text_strip_size = text_strip_size
        )
}

#' Dot plot generalized function
#' @param plotDF a data frame containing the data to be plotted.
#' @param x,y,size,fill Aesthetic mappings.
#' @param downreg_col a character string indicating the darkest color for
#' downregulated pathways. Default is `#2166AC` (blue).
#' @param upreg_col a character string indicating the darkest color for
#' upregulated pathways. Default is `#B2182B` (red).
#' @param fill_title a character string indicating the title of the fill
#' legend.
#' @param size_title a character string indicating the title of the size
#' legend.
#' @param dot_size_range a numeric vector of length 2 indicating the range of
#' dot sizes. The first value is the minimum size, and the second value is the
#' maximum size. Default is `c(0.3, 4)`.
#' @param x_position a character string indicating the position of the x-axis.
#' When `'top'`, `hjust` is set to 0. When `'bottom'`, `hjust` is set to 1.
#' @return ggplot object to be further customized.
.dotplot <- function(
        plotDF,
        x,
        y,
        size,
        fill,
        downreg_col,
        upreg_col,
        fill_title,
        size_title,
        dot_size_range,
        x_position
) {
    plotDF %>%
        ggplot(aes(
            x = x,
            y = y,
            size = size,
            fill = fill
        )) +
        geom_point(shape = 21, stroke = 0.3, color = 'black') +
        scale_fill_gradient2(
            low = downreg_col,
            mid = 'white',
            high = upreg_col,
            midpoint = 0
        ) +
        scale_size_continuous(range = dot_size_range) +
        theme_linedraw() +
        guides(
            fill = guide_colorbar(title = fill_title),
            size = guide_legend(title = size_title)
        ) +
        scale_x_discrete(position = x_position) +
        theme(
            axis.text.x = element_text(
                angle = 45,
                hjust = switch(x_position, top = 0, bottom = 1),
            ),
            axis.title.x = element_blank(),
            legend.title.position = 'left',
            legend.title = element_text(angle = 270),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank()
        )
}

#' Text size setter
#' @param text_x_size Size of cluster labels on x-axis. Default `8`.
#' @param text_y_size Size of LR-pair or pathway labels on y-axis. Default `6`.
#' @param text_strip_size Size of pathway labels in the facet strip. Default `6`.
#' @param text_title_size Size of plot title. Default `10`.
#' @param text_legend_size Size of text in legend items. Default `8`.
#' @param text_legend_title_size Size of legend title. Default `8`.
.theme_text_setter <- function(
        text_x_size = 8,
        text_y_size = 6,
        text_strip_size = 6,
        text_title_size = 10,
        text_legend_size = 8,
        text_legend_title_size = 8
) {
    theme(
        axis.text.x = element_text(size = text_x_size),
        axis.text.y = element_text(size = text_y_size),
        strip.text = element_text(size = text_strip_size),
        plot.title = element_text(size = text_title_size),
        legend.text = element_text(size = text_legend_size),
        legend.title = element_text(size = text_legend_title_size)
    )
}
