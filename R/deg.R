#' Perform pseudo-bulk differential expression test
#' @description
#' This function performs pseudo-bulk differential expression test comparing
#' two conditions (e.g. treated vs. control) within each identity class (e.g.
#' cell type or cluster). For within each identity class, we group the cells by
#' condition and biological replicate variable, and then sum the counts to
#' generate pseudo-bulk counts. The pseudo-bulk counts are then sent to DESeq2
#' for differential expression analysis.
#' @param data Full raw counts expression data of the single cells. Can be a
#' single matrix or a list of matrices without merging, a commonly seen
#' container object with raw counts at its conventional slot (e.g.
#' LayerData(seurat, 'counts'), counts(sce), or rawData(liger)).
#' @param condVar The condition variable that matches to all cells in the data.
#' @param condTest The test condition within the condition variable.
#' @param condCtrl The control condition within the condition variable.
#' @param splitVar The variable that splits the data into different identity
#' classes.
#' @param replicateVar The variable indicating the biological replicate. Default
#' `NULL`.
#' @param returnDDS Logical, whether to return the DESeqDataSet object instead
#' of just the result statistics. Default `FALSE`.
#' @param verbose Logical, whether to print some useful information. Default
#' `TRUE`
#' @param ... Arguments passed to methods and `DESeq2::DESeq()`. Note that
#' `quite = TRUE` is pre-occupied.
#' @details
#' For the input of a plain matrix or a list of matrices, any matrix coercible
#' to `dgCMatrix` is accepted. Arguments `condVar`, `splitVar` and
#' `replicateVar` must be a factor object matching all cells. For a single
#' matrix, these variables must be in the same order as the columns of the
#' matrix. For a list of matrices, the order follows as if the matrices inside
#' are concatenated in order. When presenting a list of matrices, `replicateVar`
#' is by default automatically generated following the matrix belonging.
#'
#' For container objects, `condVar`, `splitVar` and `replicateVar` are then
#' expected to be variable names in the metadata slot. For classes with
#' conventional fields in metadata, defaults are taken. For example, for a
#' Seurat object, `splitVar` is by default `Idents(data)` and `replicateVar` is
#' by default `orig.ident`. For a `liger` object, `splitVar` is by default
#' `defaultCluster` and `replicateVar` is `data$dataset`.
#' @return
#' A list of DESeq2 result objects, one for each identity class.
#' @export
#' @rdname pseudobulkDE
pseudobulkDE <- function(
        data,
        condVar,
        condTest,
        condCtrl,
        splitVar,
        replicateVar = NULL,
        ...
) {
    if (!requireNamespace('DESeq2', quietly = TRUE)) {
        cli_abort(c(
            x = 'Package {.pkg DESeq2} is required for this function.',
            i = 'Please install it with {.code BiocManager::install("DESeq2")}'
        ))
    }
    UseMethod('pseudobulkDE', data)
}

#' @export
#' @rdname pseudobulkDE
#' @method pseudobulkDE default
pseudobulkDE.default <- function(
        data,
        condVar,
        condTest,
        condCtrl,
        splitVar,
        replicateVar = NULL,
        returnDDS = FALSE,
        verbose = TRUE,
        ...
) {
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

    if (is.null(replicateVar)) {
        if (!is.list(data)) {
            cli_abort(
                'Replicate variable ({.field replicateVar}) is required when data is not a list of matrcies.'
            )
        } else {
            rep_id <- names(data) %||% as.character(seq_along(data))
            replicateVar <- rep(rep_id, sapply(data, ncol))
        }
    }
    if (length(replicateVar) != ncell) {
        cli_abort(c(
            'Length of the replicate variable ({.field replicateVar}) must be equal to the number of cells in {.field data}.',
            'i' = 'Length of replicate variable: {length(replicateVar)}',
            'i' = 'Number of columns in data matrix: {ncell}'
        ))
    }

    unique_splits <- as.character(unique(splitVar))
    result_list <- list()
    psdbulk_var <- interaction(condVar, replicateVar, drop = TRUE)
    for (i in seq_along(unique_splits)) {
        split_name <- unique_splits[i]
        if (isTRUE(verbose)) cli_process_start('Working within {.val {split_name}}')
        # figure out the number of pseudobulks to create
        cellpool <- splitVar == split_name
        # check if a condition do not contain more that one replicates

        if (inherits(data, 'dgCMatrix')) {
            psdbulk_var[!cellpool] <- NA
            psdbulk_var_design <- Matrix::fac2sparse(psdbulk_var)
            n_cells_summary <- rowSums(psdbulk_var_design)
            pseudobulk <- data %*% t(psdbulk_var_design)
            pseudobulk <- as.matrix(pseudobulk)
        } else {
            isec_genes <- Reduce(intersect, lapply(data, rownames))
            pseudobulk <- matrix(0, nrow = length(isec_genes), ncol = nlevels(psdbulk_var))
            n_cells_summary <- rep(0, nlevels(psdbulk_var))
            rep_cell_idx_start <- 0
            rep_cell_idx_end <- 0
            for (j in seq_along(data)) {
                rep_cell_idx_start <- rep_cell_idx_end + 1
                rep_cell_idx_end <- rep_cell_idx_end + ncol(data[[j]])
                rep_cell_idx <- rep_cell_idx_start:rep_cell_idx_end
                gene_fetcher <- match(isec_genes, rownames(data[[j]]))
                psdbulk_var_j <- psdbulk_var[rep_cell_idx, drop = FALSE]
                psdbulk_var_j[rep_cell_idx %in% which(!cellpool)] <- NA
                psdbulk_var_j_design <- Matrix::fac2sparse(psdbulk_var_j, drop.unused.levels = FALSE)
                n_cells_summary <- n_cells_summary + rowSums(psdbulk_var_j_design)
                pseudobulk_j <- data[[j]] %*% t(psdbulk_var_j_design)
                pseudobulk <- pseudobulk + as.matrix(pseudobulk_j)[gene_fetcher, , drop = FALSE]
            }
            dimnames(pseudobulk) <- list(isec_genes, levels(psdbulk_var))
        }
        # Make condition var on pseudobulk
        psdbulk_condVar <- tapply(
            X = as.character(condVar),
            INDEX =  psdbulk_var,
            FUN =  `[`, 1
        )
        psdbulk_condVar <- factor(
            psdbulk_condVar[colnames(pseudobulk)],
            levels = c(condCtrl, condTest)
        )
        psdbulk_repVar <- tapply(
            X = as.character(replicateVar),
            INDEX = psdbulk_var,
            FUN =  `[`, 1
        )
        meta <- data.frame(
            row.names = colnames(pseudobulk),
            condition = psdbulk_condVar,
            replicate = psdbulk_repVar[colnames(pseudobulk)],
            n_cells = n_cells_summary
        )
        cond_nrep <- meta %>%
            group_by(.data[['condition']]) %>%
            summarise(n_replicates = n_distinct(.data[['replicate']], na.rm = TRUE))
        if (any(cond_nrep$n_replicates < 2)) {
            cli_alert_danger(
                'Not all condition have at least 2 replicates. DESeq2 is very likely to fail.'
            )
            cli_alert_info('Summary of the pseudobulk generation:')
            print(meta)
        }
        result <- tryCatch(
            {
                suppressMessages({
                    dds <- DESeq2::DESeqDataSetFromMatrix(
                        countData = pseudobulk,
                        colData = meta,
                        design = stats::as.formula('~ condition')
                    )
                    dds <- DESeq2::DESeq(dds, quiet = TRUE, ...)
                })
                if (isTRUE(verbose)) cli_process_done()
                if (isTRUE(returnDDS)) dds
                else DESeq2::results(dds)
            },
            error = function(e) {
                if (isTRUE(verbose)) cli_process_failed()
                cli_alert_info('Error message:')
                cli_alert_info(e$message)
                if (isTRUE(returnDDS)) return(dds)
                else return(NULL)
            }
        )
        result_list[[split_name]] <- result
    }

    return(result_list)
}


#' @export
#' @rdname pseudobulkDE
#' @method pseudobulkDE liger
pseudobulkDE.liger <- function(
        data,
        condVar,
        condTest,
        condCtrl,
        splitVar = NULL,
        replicateVar = NULL,
        returnDDS = FALSE,
        verbose = TRUE,
        ...
) {
    if (!requireNamespace('rliger', quietly = TRUE) ||
        utils::packageVersion('rliger') < package_version('2.0.0')) {
        cli_abort(c(
            x = 'Package {.pkg rliger >= 2.0.0} is required for this function.',
            i = 'Please install it with {.code install.packages("rliger")}'
        ))
    }
    if (length(condVar) == 1) {
        condVar <- arg_match(condVar, colnames(rliger::cellMeta(data)), multiple = FALSE)
        condVar <- rliger::cellMeta(data)[[condVar]]
    }

    splitVar <- splitVar %||% rliger::defaultCluster(data)
    if (length(splitVar) == 1) {
        splitVar <- arg_match(splitVar, colnames(rliger::cellMeta(data)), multiple = FALSE)
        splitVar <- rliger::cellMeta(data)[[splitVar]]
    }

    replicateVar <- replicateVar %||% data$dataset
    if (length(replicateVar) == 1) {
        replicateVar <- arg_match(replicateVar, colnames(rliger::cellMeta(data)), multiple = FALSE)
        replicateVar <- rliger::cellMeta(data)[[replicateVar]]
    }

    data <- rliger::rawData(data)
    pseudobulkDE(
        data = data,
        condVar = condVar,
        condTest = condTest,
        condCtrl = condCtrl,
        splitVar = splitVar,
        replicateVar = replicateVar,
        returnDDS = returnDDS,
        verbose = verbose,
        ...
    )
}

#' @export
#' @rdname pseudobulkDE
#' @method pseudobulkDE Seurat
#' @param assay Name of assay to fetch counts from a Seurat object. Default
#' `"RNA"`.
pseudobulkDE.Seurat <- function(
        data,
        condVar,
        condTest,
        condCtrl,
        splitVar = NULL,
        replicateVar = NULL,
        returnDDS = FALSE,
        verbose = TRUE,
        assay = 'RNA',
        ...
) {
    if (!requireNamespace('SeuratObject', quietly = TRUE) ||
        utils::packageVersion('SeuratObject') < package_version('5.0.0')) {
        cli_abort(c(
            x = 'Package {.pkg SeuratObject >= 5.0.0} is required for this function.',
            i = 'Please install it with {.code install.packages("Seurat")}'
        ))
    }

    if (length(condVar) == 1) {
        condVar <- arg_match(condVar, colnames(data[[]]), multiple = FALSE)
        condVar <- data[[condVar, drop = TRUE]]
    }

    splitVar <- splitVar %||% SeuratObject::Idents(data)
    if (length(splitVar) == 1) {
        splitVar <- arg_match(splitVar, colnames(data[[]]), multiple = FALSE)
        splitVar <- data[[splitVar, drop = TRUE]]
    }

    replicateVar <- replicateVar %||% 'orig.ident'
    if (length(replicateVar) == 1) {
        replicateVar <- arg_match(replicateVar, colnames(data[[]]), multiple = FALSE)
        replicateVar <- data[[replicateVar, drop = TRUE]]
    }

    layerUse <- SeuratObject::Layers(data, search = 'counts', assay = assay)
    if (length(layerUse) > 1) {
        data <- lapply(layerUse, function(x) {
            SeuratObject::LayerData(data, layer = x, assay = assay)
        })
        names(data) <- layerUse
    } else {
        data <- SeuratObject::LayerData(data, layer = layerUse, assay = assay)
    }

    pseudobulkDE(
        data = data,
        condVar = condVar,
        condTest = condTest,
        condCtrl = condCtrl,
        splitVar = splitVar,
        replicateVar = replicateVar,
        returnDDS = returnDDS,
        verbose = verbose,
        ...
    )
}

######### pseudo replicate setter, too complex for now, may come back later ##
# Argument parser for creating pseudo-replicates
# @description
# Tells `pseudobulkDE()` how many pseudo-replicates to create for each
# condition.
# @param nrep The number of pseudo-replicates to create for each condition. A
# single integer applies to both conditions. A vector of two integers applies
# to test and control conditions respectively. Use `NA` within the vector to
# indicate that the condition should not be pseudo-replicated. Default `3`.
# @return A numeric vector classed as `pseudoReplicate`.
# pseudoReplicate <- function(
#         nrepTest = 3,
#         nrepCtrl = 3
# ) {
#     class(nrep) <- 'pseudoReplicate'
#     return(nrep)
# }

# Print pseudoReplicate information
# @param x An object of class `pseudoReplicate`.
# @param ... Not used.
# @method print pseudoReplicate
# print.pseudoReplicate <- function(
#         x,
#         ...
# ) {
#     print(data.frame(
#         Test = ifelse(is.na(x[1]), '(Use real)', x[1]),
#         Control = ifelse(is.na(x[2]), '(Use real)', x[2]),
#         row.names = "Number of pseudo-replicates"
#     ))
# }
