#' Differential Expression Aware Ligand-Receptor inference
#' @export
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' DEALR infers differential ligand-receptor interactions from existing
#' differential gene expression (DGE) test results. It uses Stouffer's method to
#' combine the z-scores of ligands and receptors in pair to derive combined
#' z-scores and p-values. Equation below summarizes the z-score for an LR-pair
#' where ligand (\eqn{L}) is sent from group \eqn{i} and receptor (\eqn{R}) is
#' on group \eqn{j}. We sum up the z-score of \eqn{L} tested within group
#' \eqn{i} and the z-score of \eqn{R} tested within group \eqn{j} and divide by
#' the square root of the number of z-scores (A ligand or receptor can be a
#' complex of the product of multiple genes).
#'
#' \deqn{
#' Z_{i,j,LR} = (\sum{Z_{L_i}} + \sum{Z_{R_j}})\sqrt{n}^{-1}
#' }
#'
#' The p-value is then simply the probability that the combined z-score is
#' larger than the observed value, assuming mean 0 and standard deviation 1.
#'
#' This method is assumes the DGE test to be performed with DESeq2 over
#' pseudo-bulk samples. Tests should be designed to compare to conditions
#' (e.g. treated vs control) within each group of cell (e.g. cell type).
#' Therefore, the input `degList` object is a list of DE result, each for
#' a group of cell.
#'
#' The database tested to be useful came from CellChat.
#' @section Expected input DE result:
#' This tool currently only targets at scRNAseq dataset, with the study design
#' consisting of multiple conditions each having multiple replicates. The
#' prerequisite is the pseudo-bulk DESeq2 tests comparing two conditions within
#' each group of cell. For example,
#' \enumerate{
#' \item{Optionally, subset the dataset to only the two conditions of interests
#' if there are more;}
#' \item{Split the subset by cell type and do the following step for each cell
#' type;}
#' \item{Aggregate (sum up) raw counts by replicate. This builds the
#' pseudo-bulk. This naturally yields replicate-level metadata of for the
#' conditions of interests;}
#' \item{Run DESeq2 Wald test on the pseudo-bulk data with the design on the
#' condition.}
#' }
#'
#' **Make sure to have all tests contrasting the same direction**.
#'
#' The procedure above should results in a number of DESeq2 result data frames.
#' This function expect a named list object that gathers all these data frames. List
#' names are for the cell types. **Do NOT FILTER** the DESeq2 result data frames
#' before passing to this function.
#' @param degList A named list of DESeq2 result data frames. See requirement in
#' detailed section below.
#' @param db A data frame containing ligand-receptor pair annotation. With each
#' row representing an LR-pair, it must have 1. unique rownames as for
#' identifying each LR-pair; 2. a column containing gene symbols
#' for the ligand; 3. a column containing gene symbols for the receptor; 4. a
#' column annotating the pathway the LR-pair belongs to. Complex molecules
#' transcribed from multiple genes should be presented with simple consistent
#' separator (e.g. 'Tgfbr1, Tgfbr2').
#' @param baseMeanThresh Numeric, threshold on the mean DESeq2-normalized
#' expression. Lower expression will be considered as not having sufficient
#' expression to be a ligand or receptor. Default `100`.
#' @param abslogfcThresh Numeric, threshold on the absolute value of log2 fold
#' change. Only LR-pairs with at least one ligand or receptor having large
#' enough log2 fold change will be considered. Default `1`.
#' @param baseMean_field,logfc_field,stat_field,padj_field Column names to fetch
#' necessary statistics from each DE dataframe in `degList`. Default
#' `'baseMean'`, `'log2FoldChange'`, `'stat'`, and `'padj'`, respectively.
#' @param pathway_name_field,ligand_symbol_field,receptor_symbol_field Column
#' names to fetch necessary annotation from the LR-pair database `db`.
#' Default `'pathway_name'`, `'ligand.symbol'`, and
#' `'receptor.symbol'`, respectively.
#' @param symbol_splitBy The separator applied to process the ligand and
#' receptor symbols. Default `', '`.
#' @return A `dealr` object, inherited from `tbl_df`, `dbl` and
#' `data.frame`. It contains the following columns:
#' \itemize{
#' \item{`LR_pair` - Identifier of the LR-pair, extracted from
#' `rownames(db)`}
#' \item{`sender` - Group name where the ligand statistics are from}
#' \item{`receiver` - Group name where the receptor statistics are from}
#' \item{`stat` - Combined z-score of the LR-pair between the sender and
#' receiver}
#' \item{`p` - P-value of the combined z-score}
#' \item{`pathway` - Pathway name of the LR-pair, extracted from
#' `db[[pathway_name_field]]`}
#' \item{`pathway_size` - Number of LR-pairs in the pathway, summarized
#' from `db` but not yielded from the analysis}
#' }
#' @examples
#' db_mini <- db$mouse[db$mouse$pathway_name %in% c('IL16', 'TNF'),]
#' lr <- dealr(deg_mini, db_mini)
dealr <- function(
        degList,
        db,
        baseMeanThresh = 100,
        abslogfcThresh = 1,
        baseMean_field = 'baseMean',
        logfc_field = 'log2FoldChange',
        stat_field = 'stat',
        padj_field = 'padj',
        pathway_name_field = 'pathway_name',
        ligand_symbol_field = 'ligand.symbol',
        receptor_symbol_field = 'receptor.symbol',
        symbol_splitBy = ', '
) {
    lifecycle::signal_stage('experimental', 'dealr()')
    # Input checks
    if (!is.list(degList))
        cli_abort('{.field degList} must be a list of data.frames.')
    if (is.null(names(degList))) names(degList) <- as.character(seq_along(degList))
    for (i in seq_along(degList)) {
        tryCatch(
            {
                arg_match(baseMean_field, values = colnames(degList[[i]]))
                arg_match(logfc_field, values = colnames(degList[[i]]))
                arg_match(stat_field, values = colnames(degList[[i]]))
                arg_match(padj_field, values = colnames(degList[[i]]))
                # DESeq2 original S4Vectod::DataFrame output is 3x slower in
                # subsetting tasks below
                degList[[i]] <- as.data.frame(degList[[i]])
            },
            error = function(e) {
                msg <- e$message
                cli_abort(c(
                    x = 'Error occurred when checking {.field degList[[{i}]]}',
                    i = 'Original error message: {msg}'
                ))
            }
        )
    }
    pathway_name_field <- arg_match(pathway_name_field, values = colnames(db))
    ligand_symbol_field <- arg_match(ligand_symbol_field, values = colnames(db))
    receptor_symbol_field <- arg_match(receptor_symbol_field, values = colnames(db))

    # Initialize result holder 3d-array
    # combined_z[i, j, k] fetches the combined z-score of LR-pair(k) where
    # ligands are sent from cluster(i) and receptors receives at cluster(j).
    clusters <- names(degList)
    combined_z <- array(
        data = 0,
        dim = c(length(clusters), length(clusters), nrow(db)),
        dimnames = list(clusters, clusters, rownames(db))
    )
    combined_logfc <- array(
        data = 0,
        dim = c(length(clusters), length(clusters), nrow(db)),
        dimnames = list(clusters, clusters, rownames(db))
    )

    # Go through each LR-pair, then each sender cluster for ligand(s), and then
    # each receiver cluster for receptor(s).
    cli_progress_bar(
        name = "Combining z-scores for LR-pairs",
        type = 'iterator',
        total = nrow(db) * length(clusters) * length(clusters)
    )
    for (i in seq_len(nrow(db))) {
        ligands <- strsplit(
            x = db[[ligand_symbol_field]][i],
            split = symbol_splitBy
        )[[1]]
        receptors <- strsplit(
            x = db[[receptor_symbol_field]][i],
            split =  symbol_splitBy
        )[[1]]
        for (j in seq_along(clusters)) {
            cluster_lig <- clusters[j]
            if (!all(ligands %in% rownames(degList[[cluster_lig]]))) {
                for (k in seq_along(clusters)) cli_progress_update()
                next
            }
            z_ligands <- degList[[cluster_lig]][ligands, stat_field]
            padj_ligands <- degList[[cluster_lig]][ligands, padj_field]
            baseMean_ligands <- degList[[cluster_lig]][ligands, baseMean_field]
            logfc_ligands <- degList[[cluster_lig]][ligands, logfc_field]
            if (any(is.na(padj_ligands)) ||
                any(baseMean_ligands < baseMeanThresh)) {
                for (k in seq_along(clusters)) cli_progress_update()
                next
            }
            for (k in seq_along(clusters)) {
                cluster_rec = clusters[k]
                if (!all(receptors %in% rownames(degList[[cluster_rec]]))) {
                    cli_progress_update()
                    next
                }
                z_receptors <- degList[[cluster_rec]][receptors, stat_field]
                padj_receptors <- degList[[cluster_rec]][receptors, padj_field]
                baseMean_receptors <- degList[[cluster_rec]][receptors, baseMean_field]
                logfc_receptors <- degList[[cluster_rec]][receptors, logfc_field]
                if (any(is.na(padj_receptors)) ||
                    any(baseMean_receptors < baseMeanThresh)) {
                    cli_progress_update()
                    next
                }
                if (!any(abs(c(logfc_ligands, logfc_receptors)) > abslogfcThresh)) {
                    cli_progress_update()
                    next
                }
                combined_z[j, k, i] <- stouffer_z(c(z_ligands, z_receptors))
                combined_logfc[j, k, i] <- sum(c(mean(logfc_ligands), mean(logfc_receptors)))
                cli_progress_update()
            }
        }
    }

    # Pivot the 3d-array longer so each row of the data.frame represents one
    # original entry.
    # Calculate the p-value from the combined z-score at the same time.
    combined_z_df <- .pivot_3Darray_longer(combined_z,
        values_to = 'stat',
        names_to = c('sender', 'receiver', 'LR_pair')
    ) %>%
        mutate(
            p = 2 * pnorm(abs(.data[['stat']]), lower.tail = FALSE),
            pathway = db[.data[['LR_pair']], pathway_name_field]
        )
    combined_logfc_df <- .pivot_3Darray_longer(combined_logfc,
        values_to = 'logFC',
        names_to = c('sender', 'receiver', 'LR_pair')
    )
    result <- merge(
        combined_z_df,
        combined_logfc_df,
        by = c('sender', 'receiver', 'LR_pair')
    )
    pathway_size_df <- db %>%
        group_by(.data[[pathway_name_field]]) %>%
        summarise(n = n())
    result$pathway_size <- pathway_size_df[
        match(
            result$pathway,
            pathway_size_df[[pathway_name_field]]
        ),
        'n',
        drop = TRUE
    ]
    result <- result %>%
        mutate(
            ligand_symbols = db[.data[['LR_pair']], ligand_symbol_field],
            receptor_symbols = db[.data[['LR_pair']], receptor_symbol_field],
            sender = factor(as.character(.data[['sender']]), levels = clusters),
            receiver = factor(as.character(.data[['receiver']]), levels = clusters)
        ) %>%
        select(
            .data[['LR_pair']],
            .data[['sender']],
            .data[['receiver']],
            .data[['logFC']],
            .data[['stat']],
            .data[['p']],
            .data[['pathway']],
            .data[['pathway_size']],
            .data[['ligand_symbols']],
            .data[['receptor_symbols']]
        )
    class(result) <- c('dealr', class(result))
    return(result)
}

#' Show collapsed information of DEALR result
#' @param x A `dealr` object.
#' @param trunc Number of maximum list element to show.
#' @param ... Additional arguments (not used).
#' @return NULL
#' @export
#' @method print dealr
#' @examples
#' db_mini <- db$mouse[db$mouse$pathway_name %in% c('IL16', 'TNF'),]
#' lr <- dealr(deg_mini, db_mini)
#' print(lr)
print.dealr <- function(x, trunc = 10, ...) {
    cat("DEALR result\n")
    unique_groups <- levels(x$sender)
    cat("Groups (")
    cat(length(unique_groups))
    cat("):", ansi_collapse(unique_groups, trunc = trunc, style = 'both'), '\n')
    unique_LR_pairs <- unique(x$LR_pair)
    cat("LR-pairs (")
    cat(length(unique_LR_pairs))
    cat("):", ansi_collapse(unique_LR_pairs, trunc = trunc, style = 'both'), '\n')
    cat("Pathways (")
    unique_pathways <- unique(x$pathway)
    cat(length(unique_pathways))
    cat("):", ansi_collapse(unique_pathways, trunc = trunc, style = 'both'), '\n')
    return(invisible(NULL))
}

#' Stouffer's method for combining z-scores
#' @details
#' NA input values are replaced with 0.
#' @param z Numeric vector of z-scores.
#' @return A numeric value representing the combined z-score.
#' @export
#' @examples
#' z <- c(1.5, 2.0, NA, 3.0)
#' stouffer_z(z)
stouffer_z <- function(z) {
    z[is.na(z)] <- 0
    sum(z) / sqrt(length(z))
}

#' Calculate pathway enrichment from DEALR result
#' @export
#' @description
#' Downstream of the DEALR analysis, this function further summarizes if any
#' pathway is significantly enriched between the sender and receiver. It
#' combines the z-scores of the LR-pairs in the same pathway to get the z-score
#' of the pathway. It also calculates the enrichment score of the pathway
#' between each pair of sender and receiver with
#' \deqn{\frac{(n\>significant\>LRpairs\>in\>pathway)/(size\>of\>pathway)}
#' {(n\>significant\>LRpairs)/(size\>of\>database)}\times{\sqrt{size\>of\>pathway}}}
#' @param dealr A `dealr` object. DO NOT FILTER.
#' @param p_thresh Numeric threshold on individual inference to determine if
#' the LR-pair is significant between the sender and receiver. Default
#' `0.01`.
#' @return
#' A `dealr_pe` object, inherited from `tbl_df`, `dbl` and
#' `data.frame`. It contains the following columns:
#' \itemize{
#' \item{`sender` - Group name where the ligand statistics are from}
#' \item{`receiver` - Group name where the receptor statistics are from}
#' \item{`pathway` - Pathway name}
#' \item{`stat` - Combined z-score of the pathway between the sender and
#' receiver}
#' \item{`overlap` - Number of significant LR-pairs in the
#' pathway, between the sender and receiver}
#' \item{`signif_n_intr` - Number of significant LR-pairs between the
#' sender and receiver}
#' \item{`pathway_size` - Number of LR-pairs in the pathway, summarized
#' from `db` but not yielded from the analysis}
#' \item{`p` - P-value of the combined z-score}
#' \item{`enrichment` - Enrichment score of the pathway}
#' }
#' @examples
#' db_mini <- db$mouse[db$mouse$pathway_name %in% c('IL16', 'TNF'),]
#' lr <- dealr(deg_mini, db_mini)
#' pe <- pathwayEnrich(lr)
pathwayEnrich <- function(
        dealr,
        p_thresh = 0.01
) {
    lifecycle::signal_stage('experimental', 'dealr()')
    if (!inherits(dealr, 'dealr'))
        cli_abort('{.field dealr} must be a {.cls dealr} object.')
    dealr_pe <- dealr %>%
        # Connect with database to have fold-enrichment
        group_by(.data[['sender']], .data[['receiver']]) %>%
        # How many LR-pairs are significant between one pair of sender and receiver
        mutate(signif_n_intr = sum(.data[['p']] < p_thresh)) %>%
        ungroup() %>%
        # Get pathway level stats then
        group_by(.data[['sender']], .data[['receiver']], .data[['pathway']]) %>%
        summarize(
            stat = stouffer_z(.data[['stat']]),
            # Have to get these at summarizing, individual LR-pair level information
            # is lost after this pipe
            overlap = sum(.data[['p']] < p_thresh),
            # Keep pathway level information
            signif_n_intr = unique(.data[['signif_n_intr']]),
            pathway_size = unique(.data[['pathway_size']])
        ) %>%
        mutate(
            p = 2 * pnorm(abs(.data[['stat']]), lower.tail = FALSE),
            enrichment = (.data[['overlap']] / sqrt(.data[['pathway_size']])) /
                (.data[['signif_n_intr']] / length(unique(dealr$LR_pair)))
        )
    class(dealr_pe) <- c('dealr_pe', class(dealr_pe))
    return(dealr_pe)
}

#' Show collapsed information of DEALR Pathway Enrichment result
#' @param x A `dealr_pe` object.
#' @param trunc Number of maximum list element to show.
#' @param ... Additional arguments (not used).
#' @return NULL
#' @export
#' @method print dealr_pe
#' @examples
#' db_mini <- db$mouse[db$mouse$pathway_name %in% c('IL16', 'TNF'),]
#' lr <- dealr(deg_mini, db_mini)
#' pe <- pathwayEnrich(lr)
#' print(pe)
print.dealr_pe <- function(x, trunc = 10, ...) {
    cat("DEALR Pathway Enrichment result\n")
    unique_groups <- unique(c(x$sender, x$receiver))
    cat("Groups (")
    cat(length(unique(c(x$sender, x$receiver))))
    cat("):", ansi_collapse(unique_groups, trunc = trunc, style = 'both'), '\n')
    unique_pathways <- unique(x$pathway)
    cat("Pathways (")
    cat(length(unique_pathways))
    cat("):", ansi_collapse(unique_pathways, trunc = trunc, style = 'both'), '\n')
    return(invisible(NULL))
}


#' Helper func: Look back to input stats
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' A helper function that might be removed in the future, mainly for debugging.
#' It takes the input `degList` and `db` and show the original input statistics
#' for the ligands and receptors of the LR-pairs in the pathway. Just to see if
#' we are making good decision.
#'
#' Basically what you want to see: 1. If the \code{baseMean} is high enough to
#' say the component is stably expressed, regardless of whether there is
#' significant regulation. Low \code{lfcSE} also helps indicating whether the
#' expression is stable. 2. If \code{padj} is NA, it means DESeq2 already
#' ignored it due to too-low expression. 3. Then see \code{stat} for verifying
#' if the calculation is correct.
#' @param degList The `degList` used by [dealr()].
#' @param db The `db` used by [dealr()].
#' @param pathway_use The name of the pathway to look at.
#' @param sender_use The name of the sender group.
#' @param receiver_use The name of the receiver group.
#' @param pathway_name_field,ligand_symbol_field,receptor_symbol_field Column
#' names to fetch necessary annotation from the LR-pair database `db`.
#' Default `'pathway_name'`, `'ligand.symbol'`, and
#' `'receptor.symbol'`, respectively.
#' @param symbol_splitBy The separator applied to process the ligand and
#' receptor symbols. Default `', '`.
#' @export
#' @return No return value. Print information to the console.
#' @examples
#' lookback_input(deg_mini, db$mouse, 'IL16', 'CD4_T', 'CD8_T')
lookback_input <- function(
        degList,
        db,
        pathway_use,
        sender_use,
        receiver_use,
        pathway_name_field = "pathway_name",
        ligand_symbol_field = "ligand.symbol",
        receptor_symbol_field = "receptor.symbol",
        symbol_splitBy = ", "
) {
    # Input checks
    if (!is.list(degList))
        cli_abort('{.field degList} must be a list of data.frames.')
    if (is.null(names(degList))) names(degList) <- as.character(seq_along(degList))
    sender_use <- arg_match(sender_use, values = names(degList))
    receiver_use <- arg_match(receiver_use, values = names(degList))
    pathway_name_field <- arg_match(pathway_name_field, values = colnames(db))
    ligand_symbol_field <- arg_match(ligand_symbol_field, values = colnames(db))
    receptor_symbol_field <- arg_match(receptor_symbol_field, values = colnames(db))


    db_pathway <- db %>%
        filter(.data[[pathway_name_field]] == pathway_use) %>%
        select(.data[[pathway_name_field]],
               .data[[ligand_symbol_field]],
               .data[[receptor_symbol_field]])
    cat('Pathway composition from database:\n\n')
    print(db_pathway)

    db_pathway %>%
        pull(.data[[ligand_symbol_field]]) %>%
        strsplit(symbol_splitBy) %>%
        unlist() %>%
        unique() -> ligands
    db_pathway %>%
        pull(.data[[receptor_symbol_field]]) %>%
        strsplit(symbol_splitBy) %>%
        unlist() %>%
        unique() -> receptors

    cat('\nSender DE stats:\n\n')
    print(as.data.frame(degList[[sender_use]])[ligands,])
    cat('\nReceiver DE stats:\n\n')
    print(as.data.frame(degList[[receiver_use]])[receptors,])

    return(invisible(NULL))
}
