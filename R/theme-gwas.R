#' Minimal GWAS theme for ggplot2
#'
#' A clean, publication-ready theme optimized for GWAS plots.
#'
#' @param base_size Base font size (default 11).
#' @param base_family Base font family.
#' @return A ggplot2 theme object.
#' @export
#' @examples
#' data(example_gwas)
#' manhattan_plot(example_gwas) + theme_gwas()
theme_gwas <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.y = element_line(colour = "grey30", linewidth = 0.3),
      axis.ticks.y = element_line(colour = "grey30", linewidth = 0.3),
      axis.text.x = element_text(size = base_size * 0.6, vjust = 0.5),
      plot.title = element_text(face = "bold", size = base_size * 1.2),
      legend.position = "none"
    )
}

#' Chromosome color scale
#'
#' Alternating color scale for chromosomes in Manhattan plots.
#'
#' @param colors Two-element character vector of alternating colors.
#' @param chromosomes Chromosome codes present in the data, in the order they
#'   appear along the genome. Colors alternate across this set, so any organism
#'   is handled (e.g. cattle with 29 autosomes plus X coded as 30). When `NULL`,
#'   defaults to the human coding (autosomes 1-22 plus X/Y/XY/MT as 23-26).
#' @param ... Additional arguments passed to [ggplot2::scale_color_manual()].
#' @return A ggplot2 color scale.
#' @export
#' @examples
#' data(example_gwas)
#' manhattan_plot(example_gwas, colors = c("#E64B35", "#4DBBD5"))
scale_color_chromosome <- function(colors = c("#1A5276", "#76D7C4"),
                                   chromosomes = NULL, ...) {
  chr_colors <- .chromosome_colors(colors, chromosomes)
  scale_color_manual(values = chr_colors, ...)
}

#' @rdname scale_color_chromosome
#' @export
scale_fill_chromosome <- function(colors = c("#1A5276", "#76D7C4"),
                                  chromosomes = NULL, ...) {
  chr_colors <- .chromosome_colors(colors, chromosomes)
  scale_fill_manual(values = chr_colors, ...)
}

# Build the named colour vector for a chromosome scale. Colours alternate in the
# order chromosomes are supplied (genomic order), so the number of chromosomes
# is taken from the data rather than assumed to be human (26).
.chromosome_colors <- function(colors, chromosomes) {
  if (is.null(chromosomes)) {
    chromosomes <- seq_len(26)
  }
  chromosomes <- unique(chromosomes)
  chr_colors <- rep_len(colors, length(chromosomes))
  names(chr_colors) <- as.character(chromosomes)
  chr_colors
}
