
#' Mean Functions for Gene Expression
#'
#' Functions for calculating mean gene expression across cells, designed for
#' use with [scmean()].
#'
#' @param x A count matrix with genes in rows and cells in columns.
#' @returns Numeric vector of mean values per gene.
#'
#' @details
#' `logmean` computes `rowMeans(log2(x + 1))`, the standard approach for
#' single-cell data that accounts for zero-inflation and reduces skewness.
#'
#' `trimmean` computes a trimmed mean excluding the top and bottom 5% of values
#' per row, which helps reduce the influence of outliers. Requires the `Rfast2`
#' package. When using `trimmean`, set `postFUN = log2s` in [scmean()].
#'
#' `log2s` applies the log2(x + 1) transformation.
#' @importFrom DelayedArray rowMeans
#' @export

logmean <- function(x) rowMeans(log2(x + 1))

#' @rdname logmean
#' @export
trimmean <- function(x) {
  if (!requireNamespace("Rfast2", quietly = TRUE)) {
    stop("Package 'Rfast2' is required for trimmean", call. = FALSE)
  }
  tm <- Rfast2::rowTrimMean(x)
  names(tm) <- rownames(x)
  tm
}

#' @rdname logmean
#' @export
log2s <- function(x) log2(x + 1)
