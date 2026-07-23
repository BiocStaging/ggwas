test_that("column detection is case-insensitive", {
  df <- data.frame(
    "#chrom" = 1L, pos = 1000L, rsid = "rs1", ref = "A", alt = "G",
    pvalue = 0.5, beta = 0.1, check.names = FALSE
  )
  gd <- as_gwas_data(df)
  expect_true(all(c("CHR", "BP", "P", "SNP") %in% names(gd)))
})

test_that("neg_log_pvalue / LOG10P columns are read as -log10(p)", {
  df <- data.frame(chrom = 1L, pos = 1000L, neg_log_pvalue = 8, beta = 0.1)
  gd <- as_gwas_data(df)
  expect_equal(gd$P, 1e-8)
})

test_that("read_gwas_table handles Broad-style lowercase headers", {
  f <- tempfile(fileext = ".tsv")
  writeLines(c(
    "#chrom\tpos\trsid\tref\talt\tneg_log_pvalue\tbeta\tstderr_beta\talt_allele_freq",
    "1\t727242\trs1\tG\tA\t0.66\t-0.089\t0.072\t0.14",
    "2\t758351\trs2\tA\tG\t8.5\t-0.084\t0.069\t0.15"
  ), f)
  gd <- read_gwas_table(f)
  expect_s3_class(gd, "gwas_data")
  expect_equal(nrow(gd), 2)
  expect_true(all(gd$P > 0 & gd$P <= 1))
})

test_that("missing required columns give a clear error (not a cli crash)", {
  bad <- data.frame(foo = 1:3, bar = 4:6, baz = 7:9)
  expect_error(as_gwas_data(bad), "Cannot detect required column")
})
