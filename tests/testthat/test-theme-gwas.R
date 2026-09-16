test_that("theme_gwas returns a theme", {
  expect_s3_class(theme_gwas(), "theme")
})

test_that("theme_gwas respects base_size", {
  t1 <- theme_gwas(base_size = 8)
  t2 <- theme_gwas(base_size = 14)
  expect_s3_class(t1, "theme")
  expect_s3_class(t2, "theme")
})

test_that("scale_color_chromosome returns a scale", {
  s <- scale_color_chromosome()
  expect_s3_class(s, "Scale")
})

test_that("scale_fill_chromosome returns a scale", {
  s <- scale_fill_chromosome()
  expect_s3_class(s, "Scale")
})

test_that("scale_color_chromosome accepts custom colors", {
  s <- scale_color_chromosome(colors = c("red", "blue"))
  expect_s3_class(s, "Scale")
})

test_that("chromosome scales cover non-human karyotypes", {
  # cattle: 29 autosomes plus X coded as 30. Every chromosome must get a
  # colour, including those above the old human cap of 26.
  cattle <- 1:30
  cols <- .chromosome_colors(c("#1A5276", "#76D7C4"), cattle)
  expect_equal(names(cols), as.character(cattle))
  expect_false(any(is.na(cols)))
  # colours alternate in the supplied (genomic) order
  expect_equal(unname(cols[c(1, 3, 29)]), rep("#1A5276", 3))
  expect_equal(unname(cols[c(2, 30)]), rep("#76D7C4", 2))
})

test_that("cattle Manhattan colours all chromosomes", {
  df <- do.call(rbind, lapply(1:30, function(c) {
    data.frame(SNP = paste0("s", c, "_", 1:20), CHR = c,
               BP = seq_len(20) * 1e5, P = seq(0.01, 0.99, length.out = 20))
  }))
  g <- manhattan_plot(df, genome_wide = NULL, suggestive = NULL,
                      downsample = FALSE)
  pts <- ggplot2::ggplot_build(g)$data[[1]]
  expect_false(any(is.na(pts$colour)))
})
