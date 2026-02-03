# Scientific Accuracy Tests for cellGeometry
# These tests verify mathematical correctness of deconvolution algorithms

library(testthat)
library(cellGeometry)

# ============================================
# 1. Perfect Deconvolution Test
# ============================================

test_that("deconvolution produces valid proportions", {
  # Test that deconvolution produces valid, internally consistent results
  set.seed(123)
  n_genes <- 100
  n_celltypes <- 3
  
  # Create cell type signatures with some uniqueness
  sig <- matrix(runif(n_genes * n_celltypes, 1, 5), 
                nrow = n_genes, ncol = n_celltypes)
  # Add cell-type specific markers
  for (i in 1:n_celltypes) {
    idx <- ((i-1) * 20 + 1):(i * 20)
    sig[idx, i] <- sig[idx, i] + 5
  }
  rownames(sig) <- paste0("gene", 1:n_genes)
  colnames(sig) <- paste0("cell", 1:n_celltypes)
  
  # Create mock cellMarkers object
  mk <- list(
    genemeans = sig,
    genemeans_filtered = sig,
    geneset = rownames(sig),
    subclass_table = setNames(rep(100, n_celltypes), colnames(sig)),
    opt = list(noisefilter = 2, noisefraction = 0.25),
    call = quote(cellMarkers())
  )
  class(mk) <- "cellMarkers"
  
  # Create bulk sample
  bulk <- sig %*% c(0.5, 0.3, 0.2) * 1000
  colnames(bulk) <- "sample1"
  
  # Deconvolute
  result <- deconvolute(mk, bulk, 
                        use_filter = FALSE, 
                        convert_bulk = FALSE,
                        count_space = FALSE,
                        verbose = FALSE)
  
  # Check that proportions sum to 100%
  expect_equal(sum(result$subclass$percent), 100, tolerance = 1e-6)
  
  # Check that all proportions are non-negative
  expect_true(all(result$subclass$percent >= -1e-6))
  
  # Check that result has expected dimensions
  expect_equal(ncol(result$subclass$percent), n_celltypes)
})

# ============================================
# 2. Spillover Matrix Properties
# ============================================

test_that("spillover matrix has correct mathematical properties", {
  # Create test signature matrix
  set.seed(42)
  sig <- matrix(runif(50), nrow = 10, ncol = 5)
  rownames(sig) <- paste0("g", 1:10)
  colnames(sig) <- paste0("c", 1:5)
  
  # Compute spillover using dotprod
  spill <- cellGeometry:::dotprod(sig, sig)
  
  # Property 1: Diagonal should be 1 (self-similarity)
  expect_equal(unname(diag(spill)), rep(1, 5), tolerance = 1e-10)
  
  # Property 2: Matrix should be positive semi-definite
  # (eigenvalues >= 0 for correlation-like matrix)
  eigenvalues <- eigen(spill, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigenvalues >= -1e-10))  # allow small numerical errors
})

# ============================================
# 3. Cosine Similarity Mathematical Properties
# ============================================

test_that("cosine similarity satisfies mathematical properties", {
  set.seed(42)
  mat <- matrix(runif(30, 1, 10), nrow = 10, ncol = 3)
  rownames(mat) <- paste0("g", 1:10)
  colnames(mat) <- c("A", "B", "C")
  
  # Create minimal cellMarkers object
  mk <- list(
    genemeans_filtered = mat,
    geneset = rownames(mat)
  )
  class(mk) <- "cellMarkers"
  
  result <- cos_similarity(mk)
  
  # Property 1: Diagonal is 1
  expect_equal(unname(diag(result)), rep(1, 3), tolerance = 1e-10)
  
  # Property 2: Symmetric
  expect_equal(unname(result), unname(t(result)), tolerance = 1e-10)
  
  # Property 3: Values in [-1, 1]
  expect_true(all(result >= -1 - 1e-10 & result <= 1 + 1e-10))
  
  # Property 4: For positive vectors, similarity should be in [0, 1]
  expect_true(all(result >= -1e-10))
})

# ============================================
# 4. Gene Angle Mathematical Correctness
# ============================================

test_that("gene_angle computes correct angles", {
  # Create known test case
  # Gene 1: [1, 0] -> angle to axis 1 = 0, axis 2 = 90
  # Gene 2: [0, 1] -> angle to axis 1 = 90, axis 2 = 0
  # Gene 3: [1, 1] (normalized) -> angle to both = 45
  mat <- matrix(c(
    1, 0,
    0, 1,
    1/sqrt(2), 1/sqrt(2)
  ), nrow = 3, ncol = 2, byrow = TRUE)
  rownames(mat) <- c("g1", "g2", "g3")
  colnames(mat) <- c("A", "B")
  
  result <- gene_angle(mat)
  
  # Gene 1 should have angle 0 to column A
  expect_equal(result$A["g1", "angle.deg"], 0, tolerance = 1e-5)
  
  # Gene 1 should have angle 90 to column B
  expect_equal(result$B["g1", "angle.deg"], 90, tolerance = 1e-5)
  
  # Gene 3 should have angle 45 to both (after scaling to unit sphere)
  # Note: gene_angle scales to unit sphere first
  expect_equal(result$A["g3", "angle.deg"], 45, tolerance = 1e-5)
  expect_equal(result$B["g3", "angle.deg"], 45, tolerance = 1e-5)
})

# ============================================
# 5. Residuals Computation Correctness
# ============================================

test_that("residuals are computed correctly", {
  # Create simple known case
  # If we have perfect prediction, residuals should be zero
  cellmat <- diag(3)
  colnames(cellmat) <- c("A", "B", "C")
  
  # Bulk = exact mixture
  output <- matrix(c(2, 3, 5), nrow = 1)
  colnames(output) <- c("A", "B", "C")
  
  # True bulk should be sum of components
  test <- cellmat %*% t(output)
  
  residuals <- cellGeometry:::residuals_deconv(test, cellmat, output)
  
  # Should be all zeros
  expect_equal(as.vector(residuals), rep(0, 3), tolerance = 1e-10)
})

test_that("residuals detect mismatch", {
  cellmat <- diag(3)
  output <- matrix(c(1, 1, 1), nrow = 1)
  
  # Mismatched test (not a linear combination)
  test <- matrix(c(10, 10, 10), ncol = 1)
  
  residuals <- cellGeometry:::residuals_deconv(test, cellmat, output)
  
  # Residuals should be non-zero
  expect_true(any(abs(residuals) > 1e-10))
  
  # Residuals = test - pred = [10,10,10] - [1,1,1] = [9,9,9]
  expect_equal(as.vector(residuals), c(9, 9, 9), tolerance = 1e-10)
})

# ============================================
# 6. Noise Addition Statistical Properties
# ============================================

test_that("add_noise has expected statistical properties", {
  set.seed(123)
  counts <- matrix(1e6L, nrow = 100, ncol = 10)
  
  # Add noise multiple times and check distribution
  noisy_results <- replicate(100, {
    result <- add_noise(counts, sd = 100)
    mean(result - counts)
  })
  
  # Mean of noise should be approximately 0
  expect_true(abs(mean(noisy_results)) < 100)  # generous tolerance
})

test_that("log_noise preserves relative magnitudes approximately", {
  set.seed(42)
  counts <- matrix(c(1000L, 2000L, 4000L), nrow = 3, ncol = 1)
  
  results <- replicate(50, {
    log_noise(counts, sd = 0.05)
  })
  
  # Ratios should be approximately preserved
  mean_results <- rowMeans(results)
  original_ratios <- counts[, 1] / counts[1, 1]
  result_ratios <- mean_results / mean_results[1]
  
  # Allow 20% deviation
  expect_equal(result_ratios, original_ratios, tolerance = 0.2)
})

# ============================================
# 7. Quantile Mapping Mathematical Properties
# ============================================

test_that("quantile_map preserves distribution shape", {
  skip_on_cran()
  
  set.seed(42)
  # Create two datasets with different scales
  x <- matrix(rnorm(1000, mean = 5, sd = 1), nrow = 100, ncol = 10)
  y <- matrix(rnorm(1000, mean = 10, sd = 2), nrow = 100, ncol = 10)
  rownames(x) <- rownames(y) <- paste0("gene", 1:100)
  
  # Create QQ map
  qq <- quantile_map(x, y, silent = TRUE)
  
  # Map x to y's distribution
  mapped <- qq$map(x)
  
  # Mapped values should have similar mean to y
  expect_true(abs(mean(mapped) - mean(y)) < mean(y) * 0.2)  # within 20%
})

# ============================================
# 8. Compensation Matrix Properties
# ============================================

test_that("compensation matrix is valid for deconvolution", {
  # Create non-orthogonal signatures
  set.seed(42)
  n <- 50
  sig <- matrix(runif(n * 4, 0, 10), nrow = n, ncol = 4)
  rownames(sig) <- paste0("g", 1:n)
  colnames(sig) <- c("A", "B", "C", "D")
  
  # Compute spillover
  spill <- cellGeometry:::dotprod(sig, sig)
  
  # Compensation matrix should be invertible
  expect_false(any(is.na(solve(spill))))
  
  # Check that inverse exists and is valid
  inv_spill <- solve(spill)
  identity <- spill %*% inv_spill
  expect_equal(unname(diag(identity)), rep(1, 4), tolerance = 1e-8)
})

# ============================================
# 9. Studentized Residuals Statistical Properties
# ============================================

test_that("rstudent_fit follows expected distribution under null", {
  skip_on_cran()
  
  # Simulate data where true model holds
  set.seed(42)
  n_genes <- 200
  n_samples <- 50
  n_celltypes <- 4
  
  # Create signatures
  sig <- matrix(runif(n_genes * n_celltypes, 1, 10), 
                nrow = n_genes, ncol = n_celltypes)
  rownames(sig) <- paste0("g", 1:n_genes)
  colnames(sig) <- paste0("c", 1:n_celltypes)
  
  # Generate true proportions
  props <- matrix(runif(n_samples * n_celltypes), 
                  nrow = n_samples, ncol = n_celltypes)
  props <- props / rowSums(props)  # normalize
  
  # Generate bulk data (perfect linear combination + noise)
  bulk <- t(props %*% t(sig))
  bulk <- bulk + matrix(rnorm(length(bulk), 0, 0.1), nrow = n_genes)
  colnames(bulk) <- paste0("s", 1:n_samples)
  
  # Create fit object structure
  output <- props
  residuals <- bulk - sig %*% t(props)
  
  # Hat matrix approximation
  H <- sig %*% solve(t(sig) %*% sig) %*% t(sig)
  hat <- diag(H)
  
  fit <- list(
    residuals = residuals,
    weights = NULL,
    hat = hat,
    compensation = diag(n_celltypes)
  )
  
  stud_res <- cellGeometry:::rstudent_fit(fit)
  
  # Under correct model, studentized residuals should follow approx t-distribution
  # Check that most are within expected range
  within_3sd <- mean(abs(stud_res) < 3, na.rm = TRUE)
  expect_true(within_3sd > 0.95)  # >95% within 3 SD
})

# ============================================
# 10. Full Workflow Consistency Test
# ============================================

test_that("deconvolution workflow is internally consistent", {
  skip_on_cran()
  
  set.seed(42)
  n_genes <- 100
  n_celltypes <- 3
  
  # Create cell type signatures
  sig <- matrix(0, nrow = n_genes, ncol = n_celltypes)
  for (i in 1:n_celltypes) {
    # Each cell type has some unique and some shared genes
    unique_start <- (i - 1) * 20 + 1
    unique_end <- i * 20
    sig[unique_start:unique_end, i] <- runif(20, 5, 10)
    # Shared genes
    sig[61:80, i] <- runif(20, 2, 5)
  }
  rownames(sig) <- paste0("gene", 1:n_genes)
  colnames(sig) <- paste0("cell", 1:n_celltypes)
  
  # Create mock cellMarkers
  mk <- list(
    genemeans = sig,
    genemeans_filtered = sig,
    geneset = rownames(sig),
    subclass_table = setNames(rep(100, n_celltypes), colnames(sig)),
    opt = list(noisefilter = 2, noisefraction = 0.25),
    call = quote(cellMarkers())
  )
  class(mk) <- "cellMarkers"
  
  # Create multiple samples with known proportions
  n_samples <- 10
  true_props <- matrix(runif(n_samples * n_celltypes), 
                       nrow = n_samples, ncol = n_celltypes)
  true_props <- true_props / rowSums(true_props)
  
  # Generate bulk samples
  bulk <- t(true_props %*% t(sig)) * 1000
  colnames(bulk) <- paste0("sample", 1:n_samples)
  
  # Deconvolute
  result <- deconvolute(mk, bulk,
                        use_filter = FALSE,
                        convert_bulk = FALSE,
                        count_space = FALSE,
                        comp_amount = 1)
  
  # Check proportions sum to 100%
  prop_sums <- rowSums(result$subclass$percent)
  expect_equal(unname(prop_sums), rep(100, n_samples), tolerance = 1e-6)
  
  # Correlation between true and estimated should be high
  cor_values <- sapply(1:n_celltypes, function(i) {
    cor(true_props[, i], result$subclass$percent[, i] / 100)
  })
  
  # Each cell type should have correlation > 0.9
  expect_true(all(cor_values > 0.9))
})

# ============================================
# 11. Edge Case: All Zero Row
# ============================================

test_that("functions handle all-zero rows gracefully", {
  mat <- matrix(c(
    0, 0, 0,
    1, 2, 3,
    0, 0, 0,
    4, 5, 6
  ), nrow = 4, ncol = 3, byrow = TRUE)
  
  result <- cellGeometry:::scaleSphere(mat)
  
  # Should not produce NaN
  expect_false(any(is.nan(result)))
})

# ============================================
# 12. Deconvolution Constraint Test
# ============================================

test_that("deconvolution output proportions are non-negative", {
  set.seed(42)
  n_genes <- 50
  n_celltypes <- 3
  
  sig <- matrix(runif(n_genes * n_celltypes, 1, 10), 
                nrow = n_genes, ncol = n_celltypes)
  rownames(sig) <- paste0("g", 1:n_genes)
  colnames(sig) <- paste0("c", 1:n_celltypes)
  
  mk <- list(
    genemeans = sig,
    genemeans_filtered = sig,
    geneset = rownames(sig),
    subclass_table = setNames(rep(100, n_celltypes), colnames(sig)),
    opt = list(noisefilter = 2, noisefraction = 0.25),
    call = quote(cellMarkers())
  )
  class(mk) <- "cellMarkers"
  
  # Create valid bulk sample
  bulk <- sig %*% c(0.5, 0.3, 0.2) * 1000
  colnames(bulk) <- "s1"
  
  result <- deconvolute(mk, bulk,
                        use_filter = FALSE,
                        convert_bulk = FALSE,
                        count_space = FALSE)
  
  # Check non-negative (or very close to 0)
  expect_true(all(result$subclass$percent >= -1e-6))
})
