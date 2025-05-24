#' Access miscellaneous information in the results
#' @rdname dealr_utils
#' @export
#' @description
#' - `pathways()`: Returns a character vector of all unique pathways in the result.
#' - `clusters()`: Returns a character vector of all unique clusters in the result.
#' @param x A \code{dealr} or \code{dealr_pe} object.
#' @return See Description.
pathways <- function(x) {
    if (!inherits(x, 'dealr') &&
        !inherits(x, 'dealr_pe'))
        cli_warn('The object is not a {.cls dealr} or {.cls dealr_pe} object.')
    if (!'pathway' %in% colnames(x))
        cli_abort('No column called {.field pathway}.')
    unique(x$pathway)
}

#' @rdname dealr_utils
#' @export
clusters <- function(x) {
    if (!inherits(x, 'dealr') &&
        !inherits(x, 'dealr_pe'))
        cli_warn('The object is not a {.cls dealr} or {.cls dealr_pe} object.')
    if (!'sender' %in% colnames(x) ||
        !'receiver' %in% colnames(x))
        cli_abort('Corrupted object. Both columns {.field sender} and {.field receiver} are required.')
    union(unique(x$sender), unique(x$receiver))
}

# example:
# a <- array(
#     data = 1:24,
#     dim = c(2, 3, 4),
#     dimnames = list(
#         c('a', 'b'),
#         c('c', 'd', 'e'),
#         c('f', 'g', 'h', 'i')
#     )
# )
.pivot_3Darray_longer <- function(
        x,
        values_to = 'value',
        names_to = c('dim1', 'dim2', 'dim3')
) {
    stopifnot((length(values_to) == 1))
    stopifnot((length(names_to) == 3))
    dim1names <- dimnames(x)[[1]]
    dim2names <- dimnames(x)[[2]]
    dim3names <- dimnames(x)[[3]]
    m <- dim(x)[1]
    n <- dim(x)[2]
    p <- dim(x)[3]
    data <- as.vector(x)
    dim1belongs <- rep(dim1names, times = n * p)
    dim2belongs <- rep(rep(dim2names, each = m), times = p)
    dim3belongs <- rep(dim3names, each = m * n)
    res <- data.frame(
        value = data,
        dim1 = dim1belongs,
        dim2 = dim2belongs,
        dim3 = dim3belongs
    )
    colnames(res) <- c(values_to, names_to)
    return(res)
}

# Search for top N rows in the dot plot, not the top N rows in the result table.
# like LR-pair/pathway happens to be significant in multiple cell-type-pairs,
# this only counts once.
# x - pre-filtered result table
# col_look - column name containing the LR-pair/pathway names to count for
# unique items
.search_top_n_matrow <- function(x, col_look, n) {
    x <- x %>% arrange(.data[['p']])
    seen <- list()
    i = 1
    while (length(seen) <= n) {
        if (i > nrow(x)) break

        item <- x[[col_look]][i]
        if (!item %in% names(seen)) seen[[item]] <- 1
        else seen[[item]] <- seen[[item]] + 1
        i <- i + 1
    }
    if (length(seen) > n) seen[[item]] <- NULL
    x %>% filter(.data[[col_look]] %in% names(seen))
}

.check_expr_data <- function(data) {
    checklist <- tryCatch(
        {
            if (inherits(data, 'list')) {
                # Form of a list of matrix is becoming popular in integrative
                # analysis, e.g. Seurat V5 layers, rliger datasets.
                # Merging them is expensive
                isec_genes <- Reduce(intersect, lapply(data, rownames))
                data <- lapply(data, function(x) {
                    x <- methods::as(x, 'CsparseMatrix')
                    x <- x[isec_genes, , drop = FALSE]
                    return(x)
                })
                ncell <- Reduce(`+`, sapply(data, ncol))
                list(data, ncell)
            } else {
                data <- methods::as(data, Class = 'CsparseMatrix')
                ncell <- ncol(data)
                list(data, ncell)
            }
        },
        error = function(e) {
            cli_abort('{.field data} must be a list of matrices or a matrix coercible to {.cls dgCMatrix}')
            return(list())
        }
    )
    return(checklist)
}

lognorm <- function(x) {
    x@x <- log1p(x@x / rep.int(colSums(x), diff(x@p)) * 1e6)
    return(x)
}

.categorical_col <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#FFCF00", "#aa47b9", "#e67c14",
    "#e7a2b4", "#54B0E4", "#9a5831", "#BC9DCC", "#222F75", "#1B9E77",
    "#B2DF8A", "#E3BE00", "#FF6699", "#8f3c4d", "#01e1e6", "#591cc5",
    "#A6CEE3", "#CE1261", "#8CA77B", "#5E4FA2", "#08692e", "#DCF0B9",
    "#8DD3C7", "#AAAA77"
)

.make_colors <- function(var) {
    if (is.factor(var)) {
        n <- nlevels(var)
        name <- levels(var)
    } else if (is.character(var)) {
        n <- length(unique(var))
        name <- unique(var)
    } else {
        argexp <- substitute(var)
        cli_abort('Variable {.arg {argexp}} must be a factor or character vector (categorical).')
    }
    if (n <= length(.categorical_col)) cols <- .categorical_col[seq_len(n)]
    else cols <- scales::hue_pal()(n)
    names(cols) <- name
    return(cols)
}
