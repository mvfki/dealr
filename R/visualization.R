#' Dot plot showing top differential signaling LR-pairs
#' @description
#' \itemize{
#' \item{`plotSenderLR` - Focuses on one sender on top, between any
#' receiver groups at bottom}
#' \item{`plotReceiverLR` - Focuses on one receiver at the bottom, between
#' any sender groups on the top}
#' }
#'
#' The size of dots reflects significance of the LR-pairs, and the color
#' indicates the direction of the interaction.
#' @param dealr a `dealr` object. DO NOT FILTER.
#' @param sender_use a character string indicating the sender to be plotted.
#' @param receiver_use a character string indicating the receiver to be plotted.
#' @param p_thresh a numeric value indicating the p-value threshold for
#' filtering the LR-pairs. Default is `0.01`.
#' @param pathway_use Character vector of pathways to filter the LR-pairs.
#' Applied after p-value filtering but prior to top n slicing. Default
#' `NULL`.
#' @param top_n a numeric value indicating the number of top LR-pairs to be
#' plotted, ranked by p-value. Default `NULL` means all LR-pairs, usually
#' too many.
#' @param downreg_col a character string indicating the darkest color for
#' downregulated LR-pairs. Default is `#2166AC` (blue).
#' @param upreg_col a character string indicating the darkest color for
#' upregulated LR-pairs. Default is `#B2182B` (red).
#' @param dot_size_range a numeric vector of length 2 indicating the range of
#' dot sizes. The first value is the minimum size, and the second value is the
#' maximum size. Default is `c(0.3, 4)`.
#' @return ggplot object of the dot plot.
#' @rdname plotDealr
#' @export
plotSenderLR <- function(
        dealr,
        sender_use,
        p_thresh = 0.01,
        pathway_use = NULL,
        top_n = NULL,
        downreg_col = '#2166AC',
        upreg_col = '#B2182B',
        dot_size_range = c(0.3, 4)
) {
    if (!inherits(dealr, 'dealr')) {
        cli_abort('{.field dealr} must be a {.cls dealr} object')
    }
    sender_use <- arg_match(sender_use, unique(dealr$sender))
    dealr <- dealr %>% filter(
        .data[['sender']] == sender_use,
        .data[['p']] < p_thresh)
    pathway_use <- arg_match(
        pathway_use, values = unique(dealr[['pathway']]), multiple = TRUE
    )
    if (!is.null(pathway_use))
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
        ggplot(aes(
            x = .data[['receiver']],
            y = .data[['LR_pair']],
            size = -log10(.data[['p']]),
            color = .data[['stat']]
        )) +
        geom_point() +
        theme_linedraw() +
        scale_color_gradient2(
            low = downreg_col,
            mid = 'white',
            high = upreg_col,
            midpoint = 0
        ) +
        scale_size_continuous(range = dot_size_range) +
        labs(
            title = sprintf('%s sending to others', sender_use)
        ) +
        guides(
            color = guide_colorbar(title = 'Z-score'),
            size = guide_legend(title = '-log10(p-value)')
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.title.position = 'left',
            legend.title = element_text(angle = 270)
        )
}

#' @rdname plotDealr
#' @export
plotReceiverLR <- function(
        dealr,
        receiver_use,
        p_thresh = 0.01,
        pathway_use = NULL,
        top_n = NULL,
        downreg_col = '#2166AC',
        upreg_col = '#B2182B',
        dot_size_range = c(0.3, 4)
) {
    if (!inherits(dealr, 'dealr')) {
        cli_abort('{.field dealr} must be a {.cls dealr} object')
    }
    receiver_use <- arg_match(receiver_use, unique(dealr$receiver))
    dealr <- dealr %>% filter(
        .data[['receiver']] == receiver_use,
        .data[['p']] < p_thresh
    )
    pathway_use <- arg_match(
        pathway_use, values = unique(dealr[['pathway']]), multiple = TRUE
    )
    if (!is.null(pathway_use))
        dealr <- dealr %>% filter(.data[['pathway']] %in% pathway_use)
    if (!is.null(top_n)) {
        dealr <- dealr %>%
            ungroup() %>%
            arrange(.data[['p']]) %>%
            slice_head(n = top_n)
    }
    if (nrow(dealr) == 0) {
        cli_abort('No significant LR-pairs found for the receiver {receiver_use}.')
    }
    dealr %>%
        ggplot(aes(
            x = .data[['sender']],
            y = .data[['LR_pair']],
            size = -log10(.data[['p']]),
            color = .data[['stat']]
        )) +
        geom_point() +
        theme_linedraw() +
        scale_color_gradient2(
            low = downreg_col,
            mid = 'white',
            high = upreg_col,
            midpoint = 0
        ) +
        scale_size_continuous(range = dot_size_range) +
        labs(
            title = sprintf('%s receiving from others', receiver_use)
        ) +
        guides(
            color = guide_colorbar(title = 'Z-score'),
            size = guide_legend(title = '-log10(p-value)')
        ) +
        scale_x_discrete(position = 'top') +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 0),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.title.position = 'left',
            legend.title = element_text(angle = 270)
        )
}

#' Dot plot showing pathway enrichment from differential signaling LR-pairs
#' @description
#' \itemize{
#' \item{`plotSenderPE` - Focuses on one sender on top, between any
#' receiver groups at bottom}
#' \item{`plotReceiverPE` - Focuses on one receiver at the bottom, between
#' any sender groups on the top}
#' }
#' The size of dots shows the fold enrichment of the pathways, and the color
#' indicates the significance and direction of the regulation of the pathways.
#' @param dealr_pe a `dealr_pe` object. DO NOT FILTER.
#' @param sender_use a character string indicating the sender to be plotted.
#' @param receiver_use a character string indicating the receiver to be plotted.
#' @param p_thresh a numeric value indicating the p-value threshold for
#' filtering the pathways. Default `0.01`.
#' @param top_n a numeric value indicating the number of top pathways to be
#' plotted, ranked by p-value. Default `NULL` means all pathways, usually
#' not too many.
#' @param downreg_col a character string indicating the darkest color for
#' downregulated pathways. Default is `#2166AC` (blue).
#' @param upreg_col a character string indicating the darkest color for
#' upregulated pathways. Default is `#B2182B` (red).
#' @param dot_size_range a numeric vector of length 2 indicating the range of
#' dot sizes. The first value is the minimum size, and the second value is the
#' maximum size. Default is `c(0.3, 4)`.
#' @return ggplot object of the dot plot.
#' @rdname plotDealrPE
#' @export
plotSenderPE <- function(
        dealr_pe,
        sender_use,
        p_thresh = 0.01,
        top_n = NULL,
        downreg_col = '#2166AC',
        upreg_col = '#B2182B',
        dot_size_range = c(0.3, 4)
) {
    if (!inherits(dealr_pe, 'dealr_pe')) {
        cli_abort('{.field dealr_pe} must be a {.cls dealr_pe} object')
    }
    sender_use <- arg_match(sender_use, unique(dealr_pe$sender))
    dealr_pe <- dealr_pe %>% filter(
        .data[['sender']] == sender_use,
        .data[['p']] < p_thresh
    )
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
        ggplot(aes(
            x = .data[['receiver']],
            y = .data[['pathway']],
            size = log(.data[['fold_enrichment']]),
            color = -log10(.data[['p']])*sign(.data[['stat']])
        )) +
        geom_point() +
        scale_color_gradient2(
            low = downreg_col,
            mid = 'white',
            high = upreg_col,
            midpoint = 0
        ) +
        scale_size_continuous(range = dot_size_range) +
        theme_linedraw() +
        labs(
            title = sprintf('%s sending to others', sender_use),
        ) +
        guides(
            color = guide_colorbar(title = '-log10(p) x direction'),
            size = guide_legend(title = 'log(Fold Enrichment)')
        ) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank(),
            legend.title.position = 'left',
            legend.title = element_text(angle = 270),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5)
        )
}

#' @rdname plotDealrPE
#' @export
plotReceiverPE <- function(
        dealr_pe,
        receiver_use,
        p_thresh = 0.01,
        top_n = NULL,
        downreg_col = '#2166AC',
        upreg_col = '#B2182B',
        dot_size_range = c(0.3, 4)
) {
    if (!inherits(dealr_pe, 'dealr_pe')) {
        cli_abort('{.field dealr_pe} must be a {.cls dealr_pe} object')
    }
    receiver_use <- arg_match(receiver_use, unique(dealr_pe$receiver))
    dealr_pe <- dealr_pe %>% filter(
        .data[['receiver']] == receiver_use,
        .data[['p']] < p_thresh
    )
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
        ggplot(aes(
            x = .data[['sender']],
            y = .data[['pathway']],
            size = log(.data[['fold_enrichment']]),
            color = -log10(.data[['p']])*sign(.data[['stat']])
        )) +
        geom_point() +
        scale_color_gradient2(
            low = downreg_col,
            mid = 'white',
            high = upreg_col,
            midpoint = 0
        ) +
        scale_size_continuous(range = dot_size_range) +
        theme_linedraw() +
        labs(
            title = sprintf('%s receiving from others', receiver_use),
        ) +
        guides(
            color = guide_colorbar(title = '-log10(p) x direction'),
            size = guide_legend(title = 'log(Fold Enrichment)')
        ) +
        scale_x_discrete(position = 'top') +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 0),
            axis.title.x = element_blank(),
            legend.title.position = 'left',
            legend.title = element_text(angle = 270),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5)
        )
}
