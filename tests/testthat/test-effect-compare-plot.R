test_that("effect_compare_plot returns a ggplot", {
  data(example_gwas)
  g2 <- example_gwas
  g2$BETA <- g2$BETA + stats::rnorm(nrow(g2), 0, 0.02)
  p <- effect_compare_plot(example_gwas, g2, p_threshold = 1e-3)
  expect_s3_class(p, "ggplot")
})

test_that("effect_compare_plot errors without shared variants", {
  d1 <- data.frame(SNP = paste0("a", 1:3), BETA = rnorm(3), P = runif(3))
  d2 <- data.frame(SNP = paste0("b", 1:3), BETA = rnorm(3), P = runif(3))
  expect_error(effect_compare_plot(d1, d2), "shared")
})

test_that("effect_compare_plot works without SE columns", {
  d1 <- data.frame(SNP = paste0("rs", 1:5), BETA = rnorm(5), P = runif(5))
  d2 <- data.frame(SNP = paste0("rs", 1:5), BETA = rnorm(5), P = runif(5))
  expect_s3_class(effect_compare_plot(d1, d2, show_ci = TRUE), "ggplot")
})
