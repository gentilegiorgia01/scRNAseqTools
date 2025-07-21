#' Calculate Gene Expression Summary
#'
#' For each cell, calculates the number of genes with expression >= threshold UMIs
#' and adds this information to the Seurat object metadata.
#'
#' @param seurat_obj Seurat object
#' @param threshold Numeric. Minimum UMI count threshold (default: 3)
#' @param assay Character. Name of assay to use (default: NULL uses default assay)
#'
#' @return Seurat object with added metadata column 'genes_detected'
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- calculate_expression_summary(seurat_obj, threshold = 3)
#' }
#'
#' @importFrom Seurat GetAssayData DefaultAssay
calculate_expression_summary <- function(seurat_obj, threshold = 3, assay = NULL) {
  # Input validation
  if (!methods::is(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!is.numeric(threshold) || threshold < 0) {
    stop("Threshold must be a positive numeric value")
  }

  # Check if assay is specified, otherwise use default
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }

  message("Calculating gene expression summary with threshold >= ", threshold, " UMIs...")
  message("Using assay: ", assay)

  # Get count data with explicit assay specification
  counts <- tryCatch({
    Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")
  }, error = function(e) {
    stop("Error retrieving count data: ", e$message)
  })

  # Debug information
  message("Count matrix class: ", class(counts)[1])
  message("Count matrix dimensions: ", paste(dim(counts), collapse = " x "))

  # Enhanced validation for count matrix
  if (is.null(counts)) {
    stop("Count matrix is NULL - check if the assay contains data")
  }

  if (length(dim(counts)) != 2) {
    stop("Count matrix must be 2-dimensional. Current dimensions: ", length(dim(counts)))
  }

  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("The count matrix has zero rows (", nrow(counts), ") or columns (", ncol(counts), ")")
  }

  # Convert to matrix if it's a sparse matrix (which should work with colSums)
  if (inherits(counts, "dgCMatrix") || inherits(counts, "dgTMatrix")) {
    # Sparse matrices should work with colSums, but let's be explicit
    gene_counts_per_cell <- Matrix::colSums(counts >= threshold)
  } else if (is.matrix(counts)) {
    gene_counts_per_cell <- colSums(counts >= threshold)
  } else {
    # Force conversion to matrix as last resort
    message("Converting count data to matrix format...")
    counts_matrix <- as.matrix(counts)
    gene_counts_per_cell <- colSums(counts_matrix >= threshold)
  }

  # Add to metadata
  seurat_obj$genes_detected <- gene_counts_per_cell

  # Calculate summary statistics
  summary_stats <- summary(gene_counts_per_cell)

  # Store in misc slot
  seurat_obj@misc$expression_summary <- list(
    threshold = threshold,
    assay_used = assay,
    summary_stats = summary_stats,
    total_cells = length(gene_counts_per_cell)
  )

  message("Added 'genes_detected' to cell metadata")
  message("Summary: ", paste(names(summary_stats), summary_stats, sep = "=", collapse = ", "))

  return(seurat_obj)
}
#' Plot Gene Expression Distribution
#'
#' Creates a violin plot showing the distribution of genes detected per cell
#' (genes with expression >= threshold UMIs)
#'
#' @param seurat_obj Seurat object with 'genes_detected' in metadata
#' @param title Character. Plot title (default: auto-generated)
#' @param color Character. Fill color for violin plot (default: "lightblue")
#' @param show_stats Logical. Whether to show summary statistics on plot (default: TRUE)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_expression_violin(seurat_obj)
#' print(p)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_minimal
#' @importFrom ggplot2 theme element_text annotate
plot_expression_violin <- function(seurat_obj,
                                   title = NULL,
                                   color = "lightblue",
                                   show_stats = TRUE) {

  # Input validation
  if (!methods::is(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!"genes_detected" %in% colnames(seurat_obj@meta.data)) {
    stop("'genes_detected' not found in metadata. Run calculate_expression_summary() first.")
  }

  # Get data
  genes_detected <- seurat_obj$genes_detected

  # Create data frame for plotting
  plot_data <- data.frame(
    group = "Cells",
    genes_detected = genes_detected
  )

  # Generate title if not provided
  if (is.null(title)) {
    threshold <- seurat_obj@misc$expression_summary$threshold
    if (is.null(threshold)) threshold <- "≥3"  # fallback
    title <- paste("Distribution of Genes Detected per Cell\n(Expression", threshold, "UMIs)")
  }

  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = genes_detected)) +
    ggplot2::geom_violin(fill = color, alpha = 0.7, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
    ggplot2::labs(
      title = title,
      x = "",
      y = "Number of Genes Detected"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 12, face = "bold")
    )

  # Add summary statistics if requested
  if (show_stats) {
    stats <- summary(genes_detected)
    stats_text <- paste(
      paste("Median:", round(stats["Median"], 0)),
      paste("Mean:", round(stats["Mean"], 0)),
      paste("Range:", round(stats["Min."], 0), "-", round(stats["Max."], 0)),
      sep = "\n"
    )

    # Position for text annotation (top-left)
    y_max <- max(genes_detected)
    y_pos <- y_max * 0.95

    p <- p + ggplot2::annotate(
      "text",
      x = 1,
      y = y_pos,
      label = stats_text,
      hjust = 0,
      vjust = 1,
      size = 3.5,
      color = "darkblue",
      fontface = "bold"
    )
  }

  return(p)
}


#' Alternative simplified violin plot function
#'
#' A simpler version using base R or minimal dependencies
#'
#' @param seurat_obj Seurat object with 'genes_detected' in metadata
#'
#' @return ggplot object
#' @export
plot_expression_violin_simple <- function(seurat_obj) {

  if (!methods::is(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (!"genes_detected" %in% colnames(seurat_obj@meta.data)) {
    stop("'genes_detected' not found in metadata. Run calculate_expression_summary() first.")
  }

  # Simple violin plot using Seurat's built-in function
  Seurat::VlnPlot(seurat_obj,
                  features = "genes_detected",
                  pt.size = 0) +
    ggplot2::labs(
      title = "Distribution of Genes Detected per Cell",
      x = "",
      y = "Number of Genes Detected"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
}
