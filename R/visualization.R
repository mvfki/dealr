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
#' @examples
#' db_mini <- db$mouse[db$mouse$pathway_name %in% c('IL16', 'TNF'),]
#' lr <- dealr(deg_mini, db_mini)
#' plotLRPairDot(lr)
plotLRPairDot <- function(
        dealr,
        sender_use = NULL,
        receiver_use = NULL,
        focus = NULL,
        pathway_use = NULL,
        p_thresh = 0.01,
        top_n = NULL,
        downreg_col = '#2166AC',
        downreg_name = 'Control',
        upreg_col = '#B2182B',
        upreg_name = 'Test',
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
    stopifnot(length(downreg_name) <= 1)
    stopifnot(inherits(downreg_col, 'character'))
    stopifnot(length(upreg_name) <= 1)
    stopifnot(inherits(upreg_col, 'character'))

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
            .search_top_n_matrow(col_look = 'LR_pair', n = top_n)
    }
    if (nrow(dealr) == 0) {
        cli_abort('No significant LR-pairs found for the sender {sender_use}.')
    }
    fill_title <- 'logFC'
    if (!is.null(upreg_name)) fill_title <- paste0(upreg_name, ' <- ', fill_title)
    if (!is.null(downreg_name)) fill_title <- paste0(fill_title, ' -> ', downreg_name)
    dealr %>%
        .dotplot(
            x = dealr$pair,
            y = dealr$LR_pair,
            size = -log10(dealr$p),
            # fill = dealr$stat,
            fill = dealr$logFC,
            downreg_col = downreg_col,
            upreg_col = upreg_col,
            fill_title = fill_title,
            size_title = '-log10(p-value)',
            dot_size_range = dot_size_range
        ) +
        facet_grid(
            cols = (if ('role' %in% colnames(dealr)) vars(dealr$role) else NULL),
            rows = vars(dealr$pathway),
            scales = 'free', space = 'free_y'
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
#' @param size_by Choose from `'overlap'` or `'enrichment'`. If `overlap`, dots
#' are sized by the number of significant LR-pairs annotated as in the pathway.
#' If `enrichment`, dots are sized by the log of the fold enrichment score.
#' Default `overlap`.
#' @param fe_thresh a numeric value indicating the fold enrichment threshold
#' for filtering the pathways. Default `2`.
#' @param top_n a numeric value indicating the number of top pathways to be
#' plotted, ranked by p-value. Default `NULL` means all pathways, usually
#' not too many.
#' @inheritParams .dotplot
#' @inheritParams .theme_text_setter
#' @return ggplot object of the dot plot.
#' @rdname plotDealrPE
#' @export
#' @examples
#' db_mini <- db$mouse[db$mouse$pathway_name %in% c('IL16', 'TNF'),]
#' lr <- dealr(deg_mini, db_mini)
#' pe <- pathwayEnrich(lr)
#' plotPathwayEnrichDot(pe)
plotPathwayEnrichDot <- function(
        dealr_pe,
        sender_use = NULL,
        receiver_use = NULL,
        focus = NULL,
        p_thresh = 0.01,
        size_by = c('overlap', 'enrichment'),
        fe_thresh = 2,
        top_n = NULL,
        downreg_col = '#2166AC',
        downreg_name = 'Control',
        upreg_col = '#B2182B',
        upreg_name = 'Test',
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
    size_by <- arg_match(size_by, c('overlap', 'enrichment'))

    if (size_by == 'enrichment')
        dealr_pe <- dealr_pe %>% filter(.data[['enrichment']] > fe_thresh)
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
            .search_top_n_matrow(col_look = 'pathway', n = top_n)
    }
    if (nrow(dealr_pe) == 0) {
        cli_abort('No significant LR-pairs found for the sender {sender_use}.')
    }

    dealr_pe %>%
        .dotplot(
            x = dealr_pe$pair,
            y = dealr_pe$pathway,
            size = switch(size_by,
                          overlap = dealr_pe$overlap,
                          enrichment = log(dealr_pe$enrichment)),
            fill = -log10(dealr_pe$p)*sign(dealr_pe$stat),
            downreg_col = downreg_col,
            upreg_col = upreg_col,
            fill_title = '-log10(p) x direction',
            size_title = switch(size_by,
                                overlap = 'Overlap',
                                enrichment = 'log(Enrichment)'),
            dot_size_range = dot_size_range
        ) +
        (
            if (!is.null(focus)) facet_grid(cols = vars(dealr_pe$role), space = 'free_x', scales = 'free_x')  else NULL
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
        dot_size_range
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
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            axis.title.x = element_blank(),
            legend.title.position = 'left',
            legend.title = element_text(angle = 270),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.spacing.y = unit(0, "lines"),
            panel.spacing.x = unit(0, "lines"),
            strip.text.y.right = element_text(angle = 0, color = 'black'),
            strip.text.x.top = element_text(angle = 0, color = 'black'),
            strip.background = element_rect(fill = 'grey90', color = 'black', linewidth = 0.5)
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

#' Gene expression heatmap of differential interacting LR-pair
#' @description
#' It creates heatmap focusing on a specific pair of sender and receiver cell
#' groups. The heatmap is further divided to show the contrast of between the
#' conditions where the input DEGs were derived. The columns of the heatmap
#' represent cells and each row is a gene. Rows are grouped by LR-pairs in the
#' situation when complex components are involved, and are further grouped by
#' pathways.
#' @export
#' @rdname plotLRGeneHeatmap
plotLRGeneHeatmap <- function(
        data,
        dealr,
        splitVar,
        sender_use,
        receiver_use,
        condVar,
        condTest,
        condCtrl,
        pval_thresh = 0.01,
        abslogfc_thresh = 1,
        top_n = NULL,
        pathway_use = NULL,
        splitVar_col = NULL,
        condVar_col = NULL,
        text_gene_size = 6,
        text_pathway_size = 8,
        text_condition_size = 10,
        text_annotation_size = 8,
        text_legend_size = 8,
        text_legend_title_size = 10,
        left_width = grid::unit(4, "cm"),
        middle_width = grid::unit(1, "cm"),
        right_width = grid::unit(4, "cm"),
        ...
) {
    UseMethod('plotLRGeneHeatmap', data)
}

#' @export
#' @rdname plotLRGeneHeatmap
plotLRGeneHeatmap.default <- function(
        data,
        dealr,
        splitVar,
        sender_use,
        receiver_use,
        condVar,
        condTest,
        condCtrl,
        pval_thresh = 0.01,
        abslogfc_thresh = 1,
        top_n = NULL,
        pathway_use = NULL,
        splitVar_col = NULL,
        condVar_col = NULL,
        text_gene_size = 6,
        text_pathway_size = 8,
        text_condition_size = 10,
        text_annotation_size = 8,
        text_legend_size = 8,
        text_legend_title_size = 10,
        left_width = grid::unit(4, "cm"),
        middle_width = grid::unit(1, "cm"),
        right_width = grid::unit(4, "cm"),
        ...
) {
    # All the sanity checks come first

    data_check <- .check_expr_data(data)
    data <- data_check[[1]]
    ncell <- data_check[[2]]

    if (length(condVar) != ncell) {
        cli_abort(c(
            'Length of the condition variable ({.field condVar}) must be equal to the number of cells in {.field data}.',
            'i' = 'Length of condition variable: {length(condVar)}',
            'i' = 'Number of columns in data matrix: {ncell}'
        ))
    }

    condAvail <- as.character(unique(condVar))
    condTest <- arg_match(condTest, condAvail, multiple = FALSE)
    condCtrl <- arg_match(condCtrl, condAvail, multiple = FALSE)
    condVar[!condVar %in% c(condTest, condCtrl)] <- NA
    condVar <- factor(as.character(condVar), levels = c(condCtrl, condTest))

    if (length(splitVar) != ncell) {
        cli_abort(c(
            'Length of the split variable ({.field splitVar}) must be equal to the number of cells in {.field data}.',
            'i' = 'Length of split variable: {length(splitVar)}',
            'i' = 'Number of columns in data matrix: {ncell}'
        ))
    }
    senderAvail <- intersect(
        as.character(unique(splitVar)),
        as.character(unique(dealr$sender))
    )
    if (length(senderAvail) == 0) {
        cli_abort("No common sender group found between `splitVar` and `dealr$sender`. Please check if you used consistent variable through the analysis.")
    }
    sender_use <- arg_match(sender_use, senderAvail, multiple = FALSE)

    receiverAvail <- intersect(
        as.character(unique(splitVar)),
        as.character(unique(dealr$receiver))
    )
    if (length(receiverAvail) == 0) {
        cli_abort("No common receiver group found between `splitVar` and `dealr$receiver`. Please check if you used consistent variable through the analysis.")
    }
    receiver_use <- arg_match(receiver_use, receiverAvail, multiple = FALSE)


    # Then manipulate dealr results to get information to show
    dealr <- dealr %>%
        filter(
            .data[['sender']] == sender_use,
            .data[['receiver']] == receiver_use,
            .data[['p']] < pval_thresh,
            abs(.data[['logFC']]) > abslogfc_thresh
        )
    pathway_use <- pathway_use %||% as.character(unique(dealr[['pathway']]))
    pathway_use <- arg_match(
        pathway_use, values = unique(dealr[['pathway']]), multiple = TRUE
    )
    if (!is.null(top_n)) {
        dealr <- dealr %>%
            ungroup() %>%
            .search_top_n_matrow(col_look = 'LR_pair', n = top_n)
    }
    if (nrow(dealr) == 0) {
        cli_abort('No significant LR-pairs found for the sender {sender_use}.')
    }
    lr_gene_df <- dealr %>%
        filter(.data[['pathway']] %in% pathway_use) %>%
        arrange(.data[['pathway']], .data[['p']]) %>%
        separate_rows(.data[['ligand_symbols']], sep = ",\\s*") %>%
        separate_rows(.data[['receptor_symbols']], sep = ",\\s*") %>%
        select(
            .data[['ligand_symbols']],
            .data[['receptor_symbols']],
            .data[['LR_pair']],
            .data[['pathway']],
            .data[['p']],
            .data[['logFC']]
        ) %>%
        filter(
            nchar(.data[['ligand_symbols']]) > 0,
            nchar(.data[['receptor_symbols']]) > 0
        )


    # Make the matrices for the heatmap
    idx_sender_test <- splitVar == sender_use & condVar == condTest
    idx_sender_ctrl <- splitVar == sender_use & condVar == condCtrl
    idx_receiver_test <- splitVar == receiver_use & condVar == condTest
    idx_receiver_ctrl <- splitVar == receiver_use & condVar == condCtrl
    idx_all <- which(idx_sender_test | idx_sender_ctrl |
                         idx_receiver_test | idx_receiver_ctrl)
    if (is.list(data)) {
        # Doing this because I definely want to avoid cbind'ing all matrices
        rep_cell_idx_start <- 0
        rep_cell_idx_end <- 0
        merged_data <- NULL
        for (i in seq_along(data)) {
            rep_cell_idx_start <- rep_cell_idx_end + 1
            rep_cell_idx_end <- rep_cell_idx_end + ncol(data[[i]])
            rep_cell_idx <- rep_cell_idx_start:rep_cell_idx_end
            merged_data <- cbind(merged_data, data[[i]][, rep_cell_idx %in% idx_all, drop = FALSE])
        }
        data <- merged_data
    } else {
        data <- data[, idx_all, drop = FALSE]
    }
    splitVar <- splitVar[idx_all]
    condVar <- condVar[idx_all]
    data_sender <- data[lr_gene_df$ligand_symbols, splitVar == sender_use, drop = FALSE]
    dimnames(data_sender) <- NULL
    data_receiver <- data[lr_gene_df$receptor_symbols, splitVar == receiver_use, drop = FALSE]
    dimnames(data_receiver) <- NULL

    data_sender <- as.matrix(data_sender) #%>% t() %>% scale() %>% t()
    data_receiver <- as.matrix(data_receiver) #%>% t() %>% scale() %>% t()

    # Make the heatmaps
    splitVar_col <- splitVar_col %||% .make_colors(splitVar)
    if (!sender_use %in% names(splitVar_col))
        splitVar_col[sender_use] <- .categorical_col[1]
    if (!receiver_use %in% names(splitVar_col))
        splitVar_col[receiver_use] <- .categorical_col[2]
    condVar_col <- condVar_col %||% .make_colors(condVar)
    if (!condTest %in% names(condVar_col))
        condVar_col[condTest] <- .categorical_col[3]
    if (!condCtrl %in% names(condVar_col))
        condVar_col[condCtrl] <- .categorical_col[4]

    legendTitleParam <- function(title) {
        list(
            title = title,
            title_gp = grid::gpar(fontsize = text_legend_title_size, fontface = 'bold'),
            labels_gp = grid::gpar(fontsize = text_legend_size)
        )
    }
    top_ann_df_1 <- data.frame(
        condition = condVar[splitVar == sender_use],
        role = paste0('Sender: ', splitVar[splitVar == sender_use]),
        cluster = splitVar[splitVar == sender_use]
    )
    collist1 <- list(
        condition = condVar_col[as.character(unique(top_ann_df_1$condition))],
        role = splitVar_col[as.character(unique(top_ann_df_1$cluster))]
    )
    names(collist1$role) <- paste0('Sender: ', names(collist1$role))
    heatmap_sender <- ComplexHeatmap::Heatmap(
        matrix = data_sender,
        col = circlize::colorRamp2(
            c(      0,         1,         2,         4,         6,         8),
            c('white', '#FFFFCC', '#80F2B3', '#00A0FF', '#9862B2', '#BE005F')
        ),
        name = 'Gene\nExpression',
        heatmap_legend_param = legendTitleParam('Gene\nExpression'),
        row_labels = lr_gene_df$ligand_symbols,
        row_names_gp = grid::gpar(fontsize = text_gene_size),
        row_names_side = 'left',
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
            df = top_ann_df_1[, -3],
            col = collist1,
            annotation_name_side = 'left',
            annotation_name_gp = grid::gpar(fontsize = text_annotation_size),
            annotation_legend_param = list(
                condition = legendTitleParam('Condition'),
                role = legendTitleParam('Role: cluster')
            ),
            annotation_height = grid::unit(3, "mm")
        ),
        column_split = condVar[splitVar == sender_use],
        column_title_rot = 30,
        column_title_gp = grid::gpar(fontsize = text_condition_size, fontface = 'bold'),
        cluster_column_slices = FALSE,
        row_split = lr_gene_df$pathway,
        row_title_rot = 0,
        row_title_gp = grid::gpar(fontsize = text_pathway_size),
        cluster_row_slices = FALSE,
        cluster_rows = FALSE
    )
    top_ann_df_2 <- data.frame(
        condition = condVar[splitVar == receiver_use],
        role = paste0('Receiver: ', splitVar[splitVar == receiver_use]),
        cluster = splitVar[splitVar == receiver_use]
    )
    collist2 <- list(
        condition = condVar_col[as.character(unique(top_ann_df_2$condition))],
        role = splitVar_col[as.character(unique(top_ann_df_2$cluster))]
    )
    names(collist2$role) <- paste0('Receiver: ', names(collist2$role))
    heatmap_receiver <- ComplexHeatmap::Heatmap(
        matrix = data_receiver,
        col = circlize::colorRamp2(
            c(      0,         1,         2,         4,         6,         8),
            c('white', '#FFFFCC', '#80F2B3', '#00A0FF', '#9862B2', '#BE005F')
        ),
        # col = circlize::colorRamp2(c(-2, 0, 2), c('#2166AC', 'white', '#B2182B')),
        name = 'Gene\nExpression',
        heatmap_legend_param = legendTitleParam('Gene\nExpression'),
        row_labels = lr_gene_df$receptor_symbols,
        row_names_gp = grid::gpar(fontsize = text_gene_size),
        row_names_side = 'right',
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
            df = top_ann_df_2[, -3],
            col = collist2,
            annotation_name_side = 'right',
            annotation_name_gp = grid::gpar(fontsize = text_annotation_size),
            annotation_legend_param = list(
                condition = legendTitleParam('Condition'),
                role = legendTitleParam('Role: cluster')
            ),
            annotation_height = grid::unit(3, "mm")
        ),
        column_split = condVar[splitVar == receiver_use],
        column_title_rot = 30,
        column_title_gp = grid::gpar(fontsize = text_condition_size, fontface = 'bold'),
        cluster_column_slices = FALSE,
        row_split = lr_gene_df$pathway,
        cluster_row_slices = FALSE,
        cluster_rows = FALSE
    )

    upreg_col = 'red'
    downreg_col = 'blue'
    logFC_mat <- as.matrix(lr_gene_df[, 'logFC', drop = FALSE])
    if (min(logFC_mat) > 0) {
        logFC_colmap <- circlize::colorRamp2(c(0, max(logFC_mat)), c('white', upreg_col))
    }
    if (max(logFC_mat) < 0) {
        logFC_colmap <- circlize::colorRamp2(c(min(logFC_mat), 0), c(downreg_col, 'white'))
    } else {
        logFC_colmap <- circlize::colorRamp2(
            c(min(logFC_mat), 0, max(logFC_mat)),
            c(downreg_col, 'white', upreg_col)
        )
    }
    pval_mat <- as.matrix(lr_gene_df[, 'p', drop = FALSE])
    pval_mat <- -log10(pval_mat)
    pval_colmap <- circlize::colorRamp2(
        seq(min(pval_mat), max(pval_mat), length.out = 5),
        c('#FDE725', '#5DC863', '#20908D', '#3B528B', '#440154')
    )

    logFC_middle <- ComplexHeatmap::Heatmap(
        matrix = logFC_mat,
        col = logFC_colmap,
        name = 'logFC',
        column_labels = 'logFC',
        column_names_side = 'top',
        heatmap_legend_param = legendTitleParam('logFC'),
        show_row_names = FALSE,
        cluster_rows = FALSE,
        row_split = lr_gene_df$pathway,
        row_title = FALSE,
        width = grid::unit(0.5, "cm")
    )
    pval_middle <- ComplexHeatmap::Heatmap(
        matrix = pval_mat,
        col = pval_colmap,
        name = '-log10(p)',
        column_labels = '-log10(p)',
        column_names_side = 'top',
        heatmap_legend_param = legendTitleParam('-log10(p)'),
        show_row_names = FALSE,
        cluster_rows = FALSE,
        row_split = lr_gene_df$pathway,
        row_title = FALSE,
        width = grid::unit(0.5, "cm")
    )
    ComplexHeatmap::draw(
        heatmap_sender + logFC_middle + pval_middle + heatmap_receiver,
        merge_legend = TRUE,
        heatmap_legend_side = 'right',
        annotation_legend_side = 'right',
        legend_grouping = 'original',
        gap = grid::unit(c(3, 0.5, 3), 'mm')
    )
}
