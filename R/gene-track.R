#' Gene annotation track
#'
#' Create a standalone gene annotation panel that can be composed with
#' any ggwas plot using patchwork. Genes are drawn as directional bodies
#' with strand arrows and labels; when `exon_data` is supplied, each gene is
#' rendered as an intron backbone with exon boxes and a strand arrow.
#'
#' @param gene_data A data.frame with columns: chr, start, end, gene.
#'   Optional: strand ("+"/"-").
#' @param region_chr Chromosome to display (integer or string).
#' @param region_start,region_end Region boundaries in base pairs.
#' @param exon_data Optional data.frame of exons (columns chr, start, end,
#'   gene) matching genes in `gene_data` by the `gene` column. When supplied,
#'   genes are drawn with exon structure. Read with
#'   `read_gtf(path, feature_type = "exon")`.
#' @param highlight_genes Character vector of gene names to highlight.
#' @param highlight_color Color for highlighted genes.
#' @param label_size Text size for gene labels.
#' @param track_color Default color for gene bodies.
#' @param show_strand If TRUE, draw arrows indicating gene direction.
#' @param max_genes Maximum number of genes to display. The longest genes
#'   are kept when the limit is exceeded.
#' @return A ggplot object. Compose with a Manhattan or locus plot using
#'   `patchwork::wrap_plots()`.
#' @export
#' @examples
#' genes <- data.frame(
#'   chr = c(1, 1, 1), start = c(1e6, 5e6, 8e6),
#'   end = c(2e6, 6e6, 9e6), gene = c("GeneA", "GeneB", "GeneC"),
#'   strand = c("+", "-", "+")
#' )
#'
#' # Gene-body track
#' gene_track(genes, region_chr = 1, region_start = 0, region_end = 10e6)
#'
#' # With exon structure
#' exons <- data.frame(
#'   chr = 1,
#'   start = c(1.0e6, 1.5e6, 5.0e6, 5.6e6, 8.1e6),
#'   end   = c(1.1e6, 1.7e6, 5.2e6, 5.9e6, 8.4e6),
#'   gene  = c("GeneA", "GeneA", "GeneB", "GeneB", "GeneC")
#' )
#' gene_track(genes, 1, 0, 10e6, exon_data = exons)
gene_track <- function(gene_data,
                       region_chr,
                       region_start,
                       region_end,
                       exon_data = NULL,
                       highlight_genes = NULL,
                       highlight_color = "#E74C3C",
                       label_size = 2.5,
                       track_color = "#1A5276",
                       show_strand = TRUE,
                       max_genes = 50) {

  region_chr_int <- if (is.character(region_chr)) chr_to_int(region_chr) else region_chr
  gene_data$chr_int <- if (is.character(gene_data$chr)) {
    chr_to_int(gene_data$chr)
  } else {
    as.integer(gene_data$chr)
  }

  genes <- gene_data[gene_data$chr_int == region_chr_int &
                       gene_data$end >= region_start &
                       gene_data$start <= region_end, , drop = FALSE]

  if (nrow(genes) == 0) {
    return(ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = "No genes in region", color = "grey50"))
  }

  genes$start <- pmax(genes$start, region_start)
  genes$end <- pmin(genes$end, region_end)

  if (nrow(genes) > max_genes) {
    genes$length <- genes$end - genes$start
    genes <- genes[order(-genes$length), , drop = FALSE]
    genes <- utils::head(genes, max_genes)
  }

  genes$y <- .assign_gene_tracks(genes$start, genes$end)
  genes$color <- track_color
  if (!is.null(highlight_genes)) {
    genes$color[genes$gene %in% highlight_genes] <- highlight_color
  }

  has_strand <- "strand" %in% names(genes) && show_strand
  bar_h <- 0.35
  region_span <- (region_end - region_start) / 1e6
  arrow_w <- region_span * 0.012

  if (!is.null(exon_data)) {
    plt <- .gene_track_exons(genes, exon_data, region_chr_int,
                             region_start, region_end, has_strand, bar_h)
  } else {
    plt <- .gene_track_bodies(genes, region_start, region_end,
                              has_strand, bar_h, arrow_w)
  }

  label_df <- data.frame(
    x = (genes$start + genes$end) / 2e6,
    y = genes$y + bar_h / 2 + 0.18,
    gene = genes$gene,
    stringsAsFactors = FALSE
  )

  plt +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = .data$x, y = .data$y, label = .data$gene),
      size = label_size, fontface = "italic", color = "grey20",
      direction = "x", nudge_y = 0.15,
      min.segment.length = 0, segment.color = "grey70", segment.size = 0.3,
      max.overlaps = 30, box.padding = 0.2, force = 5
    ) +
    scale_x_continuous(limits = c(region_start / 1e6, region_end / 1e6)) +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.4))) +
    labs(x = NULL, y = NULL) +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.x = element_text(size = 8),
      plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5)
    )
}

# Directional gene bodies (chevron polygons) -- default rendering.
.gene_track_bodies <- function(genes, region_start, region_end,
                               has_strand, bar_h, arrow_w) {
  arrow_polys <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
    s <- genes$start[i] / 1e6
    e <- genes$end[i] / 1e6
    yi <- genes$y[i]
    top <- yi + bar_h / 2
    bot <- yi - bar_h / 2
    clr <- genes$color[i]
    strand_i <- if (has_strand) genes$strand[i] else "."

    if (strand_i == "+") {
      tip <- min(e + arrow_w, region_end / 1e6)
      data.frame(x = c(s, e, tip, e, s), y = c(top, top, yi, bot, bot),
                 id = i, fill = clr, stringsAsFactors = FALSE)
    } else if (strand_i == "-") {
      tip <- max(s - arrow_w, region_start / 1e6)
      data.frame(x = c(e, s, tip, s, e), y = c(top, top, yi, bot, bot),
                 id = i, fill = clr, stringsAsFactors = FALSE)
    } else {
      data.frame(x = c(s, e, e, s), y = c(top, top, bot, bot),
                 id = i, fill = clr, stringsAsFactors = FALSE)
    }
  }))

  ggplot() +
    ggplot2::geom_polygon(
      data = arrow_polys,
      aes(x = .data$x, y = .data$y, group = .data$id),
      fill = arrow_polys$fill, color = NA
    )
}

# Intron backbone + exon boxes + strand arrow.
.gene_track_exons <- function(genes, exon_data, region_chr_int,
                              region_start, region_end, has_strand, bar_h) {
  ex <- exon_data
  ex$chr_int <- if (is.character(ex$chr)) chr_to_int(ex$chr) else as.integer(ex$chr)
  ex <- ex[ex$chr_int == region_chr_int &
             ex$end >= region_start & ex$start <= region_end &
             ex$gene %in% genes$gene, , drop = FALSE]
  ex$start <- pmax(ex$start, region_start)
  ex$end <- pmin(ex$end, region_end)

  ymap <- stats::setNames(genes$y, genes$gene)
  cmap <- stats::setNames(genes$color, genes$gene)
  exon_h <- bar_h * 0.9

  strand_vec <- if (has_strand) genes$strand else rep(".", nrow(genes))
  intron <- data.frame(
    x = ifelse(strand_vec == "-", genes$end, genes$start) / 1e6,
    xend = ifelse(strand_vec == "-", genes$start, genes$end) / 1e6,
    y = genes$y, color = genes$color, stringsAsFactors = FALSE
  )

  plt <- ggplot() +
    geom_segment(
      data = intron,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      color = intron$color, linewidth = 0.4,
      arrow = grid::arrow(length = grid::unit(0.045, "inches"),
                          type = "closed", ends = "last")
    )

  if (nrow(ex) > 0) {
    ex$y <- ymap[ex$gene]
    ex$fill <- cmap[ex$gene]
    plt <- plt + geom_rect(
      data = ex,
      aes(xmin = .data$start / 1e6, xmax = .data$end / 1e6,
          ymin = .data$y - exon_h / 2, ymax = .data$y + exon_h / 2),
      fill = ex$fill, color = NA
    )
  }

  plt
}
