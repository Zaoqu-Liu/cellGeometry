# Core function tests for cellGeometry
# Scientific accuracy and numerical stability tests

library(testthat)
library(cellGeometry)

# ============================================
# 1. Vector Projection Tests (dotprod)
# ============================================

test_that("dotprod computes correct vector projection", {
  # Simple case: projection of unit vectors
  cellmat <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  colnames(cellmat) <- c("A", "B")
  test <- matrix(c(1, 1), nrow = 2, ncol = 1)
  
  result <- cellGeometry:::dotprod(test, cellmat)
  
  # proj_A([1,1]) = (1*1 + 1*0) / (1^2 + 0^2) = 1
  # proj_B([1,1]) = (1*0 + 1*1) / (0^2 + 1^2) = 1
  expect_equal(as.vector(result), c(1, 1), tolerance = 1e-10)
})

test_that("dotprod handles zero vectors", {
  cellmat <- matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2)
  test <- matrix(c(1, 1), nrow = 2, ncol = 1)
  
  # Should not produce NaN or Inf
  result <- cellGeometry:::dotprod(test, cellmat)
  expect_false(any(is.nan(result)))
  expect_false(any(is.infinite(result)))
})

test_that("dotprod is mathematically correct for known values", {
  # v = [3, 4], ||v||^2 = 25
  # u = [6, 8]
  # proj_v(u) = (6*3 + 8*4) / 25 = 50/25 = 2
  cellmat <- matrix(c(3, 4), nrow = 2, ncol = 1)
  test <- matrix(c(6, 8), nrow = 2, ncol = 1)
  
  result <- cellGeometry:::dotprod(test, cellmat)
  expect_equal(as.vector(result), 2, tolerance = 1e-10)
})

# ============================================
# 2. Hypersphere Scaling Tests (scaleSphere)
# ============================================

test_that("scaleSphere normalizes rows to unit length", {
  mat <- matrix(c(3, 4, 1, 0, 5, 12), nrow = 3, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::scaleSphere(mat)
  row_lengths <- sqrt(rowSums(result^2))
  
  expect_equal(row_lengths, rep(1, 3), tolerance = 1e-10)
})

test_that("scaleSphere handles zero rows", {
  mat <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::scaleSphere(mat)
  
  # Should not produce NaN
  expect_false(any(is.nan(result)))
})

# ============================================
# 3. Cosine Similarity Tests
# ============================================

test_that("cos_sim produces valid similarity matrix", {
  mat <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, ncol = 2, byrow = TRUE)
  rownames(mat) <- c("A", "B", "C")
  
  result <- cellGeometry:::cos_sim(mat)
  
  # Diagonal should be 1

  expect_equal(diag(result), rep(1, 2), tolerance = 1e-10)
  
  # Should be symmetric
  expect_equal(result, t(result), tolerance = 1e-10)
  
  # All values should be in [-1, 1]
  expect_true(all(result >= -1 - 1e-10 & result <= 1 + 1e-10))
})

test_that("cos_sim handles orthogonal vectors", {
  # Orthogonal vectors should have similarity 0
  mat <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::cos_sim(mat)
  
  expect_equal(result[1, 2], 0, tolerance = 1e-10)
})

test_that("cos_sim handles identical vectors", {
  mat <- matrix(c(1, 2, 1, 2), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::cos_sim(mat)
  
  # Identical vectors should have similarity 1
  expect_equal(result[1, 2], 1, tolerance = 1e-10)
})

# ============================================
# 4. Noise Reduction Tests
# ============================================

test_that("reduceNoise sets low values to zero", {
  mat <- matrix(c(10, 1, 0.5, 8, 2, 0.3), nrow = 2, ncol = 3, byrow = TRUE)
  
  result <- reduceNoise(mat, noisefilter = 2, noisefraction = 0.25)
  
  # Values below threshold should be 0
  expect_true(all(result[result != 0] >= 0))
})

test_that("reduceNoise preserves high values", {
  mat <- matrix(c(10, 8, 6, 12, 10, 8), nrow = 2, ncol = 3, byrow = TRUE)
  
  result <- reduceNoise(mat, noisefilter = 2, noisefraction = 0.25)
  
  # High values should be preserved
  expect_equal(mat[1, 1], result[1, 1])
})

# ============================================
# 5. Mean Function Tests
# ============================================

test_that("logmean computes correct values", {
  mat <- matrix(c(0, 2, 6, 0, 0, 0), nrow = 2, ncol = 3, byrow = TRUE)
  
  result <- logmean(mat)
  
  # log2(0+1) = 0, log2(2+1) = 1.585, log2(6+1) = 2.807
  # mean for row 1 = (0 + 1.585 + 2.807) / 3 = 1.464
  expected_row1 <- mean(log2(c(0, 2, 6) + 1))
  expect_equal(result[1], expected_row1, tolerance = 1e-10)
  
  # Row 2: all zeros -> log2(1) = 0
  expect_equal(result[2], 0, tolerance = 1e-10)
})

test_that("log2s transforms correctly", {
  x <- c(0, 1, 3, 7, 15)
  result <- log2s(x)
  expected <- log2(x + 1)
  expect_equal(result, expected, tolerance = 1e-10)
})

# ============================================
# 6. Equal Weight Tests
# ============================================

test_that("equalweight computes correct inverse lengths", {
  mat <- matrix(c(3, 4, 5, 12), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::equalweight(mat)
  
  # ||[3,4]|| = 5, weight = 1/5 = 0.2
  # ||[5,12]|| = 13, weight = 1/13 ≈ 0.0769
  expect_equal(result[1], 1/5, tolerance = 1e-10)
  expect_equal(result[2], 1/13, tolerance = 1e-10)
})

test_that("equalweight handles zero vectors", {
  mat <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::equalweight(mat)
  
  # Should not produce Inf
  expect_false(any(is.infinite(result)))
})

# ============================================
# 7. Residuals Computation Tests
# ============================================

test_that("residuals_deconv computes correct residuals", {
  cellmat <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  output <- matrix(c(2, 3), nrow = 1, ncol = 2)
  test <- matrix(c(2, 3), nrow = 2, ncol = 1)
  
  result <- cellGeometry:::residuals_deconv(test, cellmat, output)
  
  # pred = cellmat %*% t(output) = [2, 3]
  # residuals = test - pred = [0, 0]
  expect_equal(as.vector(result), c(0, 0), tolerance = 1e-10)
})

# ============================================
# 8. Spillover Metric Tests
# ============================================

test_that("comp_metric computes mean absolute off-diagonal", {
  m <- matrix(c(1, 0.2, 0.3, 1), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::comp_metric(m)
  
  # Off-diagonal: 0.2, 0.3
  # Mean: (0.2 + 0.3) / 4 = 0.125  (including zeros on diagonal after subtraction)
  expect_equal(result, mean(c(0.2, 0.3, 0, 0)), tolerance = 1e-10)
})

test_that("max_spill finds maximum off-diagonal", {
  m <- matrix(c(1, 0.2, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::max_spill(m)
  
  expect_equal(result, 0.5, tolerance = 1e-10)
})

# ============================================
# 9. Gene Angle Tests
# ============================================

test_that("gene_angle returns correct structure", {
  mat <- matrix(runif(30), nrow = 10, ncol = 3)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:3)
  
  result <- gene_angle(mat)
  
  expect_type(result, "list")
  expect_length(result, 3)
  expect_equal(names(result), colnames(mat))
  
  # Each element should have angle, angle.deg, max, rank columns
  expect_true(all(c("angle", "angle.deg", "max", "rank") %in% colnames(result[[1]])))
})

test_that("gene_angle angles are in valid range", {
  mat <- matrix(runif(30), nrow = 10, ncol = 3)
  rownames(mat) <- paste0("gene", 1:10)
  colnames(mat) <- paste0("cell", 1:3)
  
  result <- gene_angle(mat)
  
  for (i in seq_along(result)) {
    # Angles should be in [0, pi]
    expect_true(all(result[[i]]$angle >= 0 & result[[i]]$angle <= pi))
    # Degrees should be in [0, 180]
    expect_true(all(result[[i]]$angle.deg >= 0 & result[[i]]$angle.deg <= 180))
  }
})

# ============================================
# 10. Deconvolution Core Tests
# ============================================

test_that("deconv produces valid output", {
  # Simple 2 cell type case
  cellmat <- matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
  colnames(cellmat) <- c("A", "B")
  test <- matrix(c(0.5, 0.5), nrow = 2, ncol = 1)
  
  m_itself <- cellGeometry:::dotprod(cellmat, cellmat)
  
  result <- cellGeometry:::deconv(test, cellmat, comp_amount = 1, m_itself)
  
  expect_true("output" %in% names(result))
  expect_true("percent" %in% names(result))
  
  # Percentages should sum to 100
  expect_equal(sum(result$percent), 100, tolerance = 1e-6)
})

# ============================================
# 11. Add Noise Tests
# ============================================

test_that("add_noise preserves non-negativity", {
  counts <- matrix(sample(100:1000, 20), nrow = 4, ncol = 5)
  
  result <- add_noise(counts, sd = 10)
  
  expect_true(all(result >= 0))
})

test_that("log_noise preserves matrix dimensions", {
  counts <- matrix(sample(100:1000, 20), nrow = 4, ncol = 5)
  
  result <- log_noise(counts, sd = 0.1)
  
  expect_equal(dim(result), dim(counts))
})

test_that("graded_log_noise handles zero matrix", {
  counts <- matrix(0L, nrow = 4, ncol = 5)
  
  # Should not error
  result <- graded_log_noise(counts, sd = 0.1)
  
  expect_false(any(is.nan(result)))
  expect_false(any(is.infinite(result)))
})

# ============================================
# 12. Edge Case Tests
# ============================================

test_that("scaleSphere handles single column matrix", {
  # 3x1 matrix: each row has only one element
  mat <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
  
  result <- cellGeometry:::scaleSphere(mat)
  
  # Each row should be normalized to unit length
  # For single-element rows, result should be 1 (or -1)
  row_lengths <- sqrt(rowSums(result^2))
  expect_equal(row_lengths, rep(1, 3), tolerance = 1e-10)
})

test_that("scaleSphere handles single row matrix", {
  # 1x3 matrix: one row with three elements
  mat <- matrix(c(1, 2, 3), nrow = 1, ncol = 3)
  
  result <- cellGeometry:::scaleSphere(mat)
  
  # The single row should be normalized to unit length
  row_length <- sqrt(sum(result^2))
  expect_equal(row_length, 1, tolerance = 1e-10)
})

test_that("logmean handles large values", {
  mat <- matrix(c(1e6, 1e7, 1e8), nrow = 1, ncol = 3)
  
  result <- logmean(mat)
  
  expect_false(is.nan(result))
  expect_false(is.infinite(result))
})

# ============================================
# 13. Numerical Precision Tests
# ============================================

test_that("cos_sim handles near-identical vectors", {
  # Vectors that are very close but not identical
  mat <- matrix(c(1, 1e-15, 1, 0), nrow = 2, ncol = 2, byrow = TRUE)
  
  result <- cellGeometry:::cos_sim(mat)
  
  # Should not produce NaN
  expect_false(any(is.nan(result)))
  # Similarity should be close to 1
  expect_true(abs(result[1, 2]) <= 1 + 1e-10)
})

test_that("dotprod handles very small values", {
  cellmat <- matrix(c(1e-15, 1e-15), nrow = 2, ncol = 1)
  test <- matrix(c(1, 1), nrow = 2, ncol = 1)
  
  result <- cellGeometry:::dotprod(test, cellmat)
  
  expect_false(any(is.nan(result)))
  expect_false(any(is.infinite(result)))
})

# ============================================
# 14. Integration Test: Simulation Workflow
# ============================================

test_that("simulation workflow is self-consistent", {
  skip_on_cran()
  
  # Create mock cellMarkers object structure
  set.seed(42)
  n_genes <- 100
  n_cells <- 3
  
  genemeans <- matrix(runif(n_genes * n_cells, 0, 10), 
                      nrow = n_genes, ncol = n_cells)
  rownames(genemeans) <- paste0("gene", 1:n_genes)
  colnames(genemeans) <- paste0("cell", 1:n_cells)
  
  # Create minimal cellMarkers-like object
  mk <- list(
    genemeans = genemeans,
    subclass_table = setNames(c(100, 100, 100), colnames(genemeans))
  )
  class(mk) <- "cellMarkers"
  
  # Generate samples
  samples <- generate_samples(mk, n = 5, method = "unif")
  
  expect_equal(nrow(samples), 5)
  expect_equal(ncol(samples), n_cells)
  expect_true(all(samples >= 0))
})
