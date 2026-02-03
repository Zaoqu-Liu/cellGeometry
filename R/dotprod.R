
# Vector projection: proj_v(u) = (u·v) / ||v||^2
dotprod <- function(test, cellmat) {
  md <- colSums(cellmat^2)
  md[md == 0] <- .Machine$double.eps
  t(crossprod(cellmat, test) / md)
}

# Compute weights for equal gene contribution
equalweight <- function(cellmat) {

  vec_length <- sqrt(rowSums(cellmat^2))
  vec_length[vec_length == 0] <- .Machine$double.eps
  1 / vec_length
}

# Mean absolute off-diagonal spillover
comp_metric <- function(m) {
  diag(m) <- 0

  mean(abs(m))
}

# Maximum off-diagonal spillover
max_spill <- function(m) {
  diag(m) <- 0
  max(m)
}

# Extreme absolute value
max_abs <- function(m) {
  if (abs(min(m)) > max(m)) return(min(m))
  max(m)
}

# Compute residuals: observed - predicted
residuals_deconv <- function(test, cellmat, output) {
  pred <- tcrossprod(cellmat, output)
  test - pred
}
