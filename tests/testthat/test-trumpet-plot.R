test_that("trumpet_plot returns a ggplot", {
  data(example_gwas)
  p <- trumpet_plot(example_gwas, n = 50000)
  expect_s3_class(p, "ggplot")
})

test_that("trumpet_plot supports signed effects and labels", {
  data(example_gwas)
  p <- trumpet_plot(example_gwas, n = 50000, signed = TRUE, label_top_n = 3)
  expect_s3_class(p, "ggplot")
})

test_that("trumpet_plot errors when no sample size is available", {
  df <- data.frame(BETA = c(0.1, 0.2), AF = c(0.2, 0.3), P = c(1e-9, 0.5))
  expect_error(trumpet_plot(df), "sample size|n_col|n")
})

test_that("trumpet_plot power contours flare toward rare variants", {
  df <- data.frame(BETA = rnorm(200, 0, 0.1),
                   AF = runif(200, 0.01, 0.5),
                   P = runif(200))
  p <- trumpet_plot(df, n = 10000)
  expect_s3_class(p, "ggplot")
})
