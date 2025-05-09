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
