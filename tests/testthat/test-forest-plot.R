test_that("forest_plot returns a ggplot", {
  df <- data.frame(SNP = paste0("rs", 1:4),
                   BETA = c(0.2, -0.1, 0.35, 0.05),
                   SE = c(0.05, 0.04, 0.08, 0.03))
  expect_s3_class(forest_plot(df), "ggplot")
})

test_that("forest_plot exponentiate and ordering work", {
  df <- data.frame(SNP = paste0("rs", 1:4),
                   BETA = c(0.2, -0.1, 0.35, 0.05),
                   SE = c(0.05, 0.04, 0.08, 0.03))
  expect_s3_class(forest_plot(df, exponentiate = TRUE, order_by = "effect"),
                  "ggplot")
})

test_that("forest_plot supports grouping", {
  df <- data.frame(SNP = rep(paste0("rs", 1:3), 2),
                   BETA = rnorm(6), SE = runif(6, 0.02, 0.06),
                   cohort = rep(c("A", "B"), each = 3))
  expect_s3_class(forest_plot(df, group = "cohort"), "ggplot")
})

test_that("forest_plot errors on missing columns", {
  expect_error(forest_plot(data.frame(x = 1), effect = "BETA"), "not found")
})
