#' Compare variant effects between two GWAS
#'
#' Scatter the per-variant effect sizes of two studies against each other for
#' the variants they share, to assess replication and effect concordance
#' (for example discovery versus replication, or two ancestries). A `y = x`
#' reference line is drawn; points are colored by whether the two effects
#' agree in sign among variants significant in either study.
#'
#' @param gwas1,gwas2 `gwas_data` objects or data.frames. Variants are matched
#'   on the SNP column.
#' @param snp,beta,se,p Column name overrides applied to both inputs when they
#'   are plain data.frames.
#' @param labels Length-2 character vector of axis labels.
#' @param p_threshold Significance threshold used to classify variants.
#' @param show_ci If TRUE, draw confidence-interval crosses for the
#'   significant variants (requires SE columns).
#' @param ci Confidence level for the crosses.
#' @param label_top_n Label the N most significant shared variants.
#' @param colors Named colors for "concordant", "discordant" and "ns".
#' @param point_size Point size.
#' @param title Plot title.
#' @return A ggplot object.
#' @export
#' @examples
#' data(example_gwas)
#' g2 <- example_gwas
#' g2$BETA <- g2$BETA + stats::rnorm(nrow(g2), 0, 0.02)
#' effect_compare_plot(example_gwas, g2,
#'                     labels = c("Discovery", "Replication"),
#'                     p_threshold = 1e-3)
effect_compare_plot <- function(gwas1, gwas2,
                                snp = NULL, beta = NULL, se = NULL, p = NULL,
                                labels = c("Study 1", "Study 2"),
                                p_threshold = 5e-8,
                                show_ci = TRUE,
                                ci = 0.95,
                                label_top_n = 10,
                                colors = c(concordant = "#2C7FB8",
                                           discordant = "#D7301F",
                                           ns = "#BDC3C7"),
                                point_size = 1.8,
                                title = NULL) {

  prep <- function(g) {
    g <- as.data.frame(g)
    scol <- if (!is.null(snp)) snp else if ("SNP" %in% names(g)) "SNP" else NULL
    bcol <- if (!is.null(beta)) beta else if ("BETA" %in% names(g)) "BETA" else NULL
    if (is.null(scol) || is.null(bcol)) {
      cli_abort("Both datasets need SNP and BETA columns (or specify via {.arg snp} and {.arg beta}).")
    }
    secol <- if (!is.null(se)) se else if ("SE" %in% names(g)) "SE" else NULL
    pcol <- if (!is.null(p)) p else if ("P" %in% names(g)) "P" else NULL
    data.frame(
      SNP = as.character(g[[scol]]),
      b = g[[bcol]],
      se = if (!is.null(secol)) g[[secol]] else NA_real_,
      p = if (!is.null(pcol)) g[[pcol]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  d1 <- prep(gwas1)
  d2 <- prep(gwas2)
  m <- merge(d1, d2, by = "SNP", suffixes = c("1", "2"))
  m <- m[!is.na(m$b1) & !is.na(m$b2), , drop = FALSE]
  if (nrow(m) == 0) cli_abort("No shared variants matched on {.field SNP}.")

  sig <- (!is.na(m$p1) & m$p1 < p_threshold) |
         (!is.na(m$p2) & m$p2 < p_threshold)
  m$class <- ifelse(!sig, "ns",
                    ifelse(sign(m$b1) == sign(m$b2), "concordant", "discordant"))
  m$class <- factor(m$class, levels = c("concordant", "discordant", "ns"))

  plt <- ggplot(m, aes(x = .data$b1, y = .data$b2)) +
    geom_hline(yintercept = 0, color = "grey85") +
    ggplot2::geom_vline(xintercept = 0, color = "grey85") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50")

  if (show_ci && any(!is.na(m$se1)) && any(!is.na(m$se2))) {
    z <- stats::qnorm(1 - (1 - ci) / 2)
    ci_df <- m[sig & !is.na(m$se1) & !is.na(m$se2), , drop = FALSE]
    if (nrow(ci_df) > 0) {
      plt <- plt +
        geom_segment(data = ci_df,
          aes(x = .data$b1 - z * .data$se1, xend = .data$b1 + z * .data$se1,
              y = .data$b2, yend = .data$b2),
          color = "grey70", linewidth = 0.3, inherit.aes = FALSE) +
        geom_segment(data = ci_df,
          aes(x = .data$b1, xend = .data$b1,
              y = .data$b2 - z * .data$se2, yend = .data$b2 + z * .data$se2),
          color = "grey70", linewidth = 0.3, inherit.aes = FALSE)
    }
  }

  plt <- plt +
    geom_point(aes(color = .data$class), size = point_size, alpha = 0.75) +
    scale_color_manual(values = colors, name = NULL, drop = FALSE) +
    labs(x = paste(labels[1], "effect"), y = paste(labels[2], "effect"),
         title = title) +
    theme_gwas()

  r <- suppressWarnings(stats::cor(m$b1, m$b2, use = "complete.obs"))
  if (is.finite(r)) {
    plt <- plt + labs(subtitle = sprintf("Pearson r = %.2f (%d shared variants)",
                                          r, nrow(m)))
  }

  if (!is.null(label_top_n)) {
    m$minp <- pmin(m$p1, m$p2, na.rm = TRUE)
    top <- utils::head(m[order(m$minp), ], label_top_n)
    plt <- plt + ggrepel::geom_text_repel(
      data = top, aes(x = .data$b1, y = .data$b2, label = .data$SNP),
      inherit.aes = FALSE, size = 2.8, max.overlaps = 15,
      segment.color = "grey50"
    )
  }

  plt
}
