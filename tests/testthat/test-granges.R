test_that("as_granges converts gwas_data to GRanges", {
  data(example_gwas)
  gr <- as_granges(example_gwas)
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), nrow(example_gwas))
  expect_true("P" %in% names(S4Vectors::mcols(gr)))
})

test_that("as_gwas_data accepts GRanges input (round-trip)", {
  data(example_gwas)
  gr <- as_granges(example_gwas)
  gd <- as_gwas_data(gr)
  expect_s3_class(gd, "gwas_data")
  expect_equal(nrow(gd), nrow(example_gwas))
  expect_true(all(c("CHR", "BP", "P") %in% names(gd)))
})

test_that("plots accept GRanges via as_gwas_data", {
  data(example_gwas)
  gr <- as_granges(example_gwas)
  expect_s3_class(manhattan_plot(as_gwas_data(gr)), "ggplot")
})
