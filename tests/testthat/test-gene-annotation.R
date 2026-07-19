test_that("gene_annotation returns expected structure for both builds", {
  for (b in c("GRCh38", "GRCh37")) {
    g <- gene_annotation(b)
    expect_s3_class(g, "data.frame")
    expect_true(all(c("chr", "start", "end", "gene", "strand") %in% names(g)))
    expect_gt(nrow(g), 15000)
    expect_type(g$chr, "integer")
  }
})

test_that("gene_annotation rejects unknown builds", {
  expect_error(gene_annotation("hg19"))
})

test_that("built-in genes feed gene_track", {
  g <- gene_annotation("GRCh38")
  p <- gene_track(g, region_chr = 6, region_start = 25e6, region_end = 34e6,
                  max_genes = 20)
  expect_s3_class(p, "ggplot")
})
