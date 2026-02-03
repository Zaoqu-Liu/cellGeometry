
#' Reduce noise in single-cell data
#' 
#' Filter for reducing noise in single-cell mean expression data. For each gene,
#' expression values below a threshold are set to zero. The threshold is the 
#' minimum of `noisefilter` and `noisefraction * max_expression`.
#' 
#' @param cellmat Matrix of log2 mean gene expression with genes in rows and 
#'   cell types in columns.
#' @param noisefilter Numeric. Upper bound for the noise threshold. Expression
#'   values above this level are always retained.
#' @param noisefraction Numeric (0-1). Fraction of maximum expression per gene
#'   used as threshold. Values below `max_expression * noisefraction` are set
#'   to zero, unless this exceeds `noisefilter`.
#' @returns Matrix with low expression values set to zero.
#' @details
#' For each gene (row), the threshold is: `min(max_expr * noisefraction, noisefilter)`.
#' This means highly expressed genes have stricter filtering, but the threshold

#' never exceeds `noisefilter`.
#' @export
#'
reduceNoise <- function(cellmat,
                        noisefilter = 2,
                        noisefraction = 0.25) {
  genemax <- rowMaxs(cellmat)
  threshold <- pmin(genemax * noisefraction, noisefilter)
  # Compare each row element with its corresponding threshold
  cellmat[cellmat < threshold] <- 0
  cellmat
}
