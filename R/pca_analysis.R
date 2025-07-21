#' Perform PCA Analysis and Plot Variance Explained
#'
#' Normalizes, scales and performs PCA on a Seurat object. Saves PCA variance plot.
#'
#' @param seurat_obj A Seurat object
#' @param n_pcs Integer. Number of principal components to compute (default: 20)
#' @param save_path_rds Path to save the Seurat object after PCA (default: NULL)
#' @param save_plot Path to save PCA variance histogram as PNG (default: NULL)
#'
#' @return A Seurat object with PCA computed
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- perform_pca_analysis(
#'   seurat_obj,
#'   n_pcs = 20,
#'   save_path_rds = "data/seurat_obj.pca.rds",
#'   save_plot = "data/pca_variance_histogram.png"
#' )
#' }
#'
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA VariableFeatures
#' @importFrom ggplot2 ggplot aes geom_bar labs theme_minimal scale_x_continuous
perform_pca_analysis <- function(seurat_obj,
                                 n_pcs = 20,
                                 save_path_rds = NULL,
                                 save_plot = NULL) {

  if (!methods::is(seurat_obj, "Seurat")) stop("Input must be a Seurat object")
  if (!is.numeric(n_pcs) || n_pcs <= 0) stop("n_pcs must be a positive integer")

  message("Performing PCA on the dataset...")

  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 2000)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(seurat_obj), npcs = n_pcs)

  # Variance explained
  pca_variance <- (seurat_obj@reductions$pca@stdev)^2
  pca_variance_prop <- pca_variance / sum(pca_variance)

  pca_df <- data.frame(
    PC = 1:n_pcs,
    Variance_Explained = pca_variance_prop[1:n_pcs]
  )

  # Plot histogram
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC, y = Variance_Explained)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    ggplot2::labs(title = "Variance Explained by Principal Components",
                  x = "Principal Component",
                  y = "Proportion of Variance Explained") +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = 1:n_pcs)

  # Save plot if path provided
  if (!is.null(save_plot)) {
    dir.create(dirname(save_plot), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_plot, plot = p, width = 6, height = 4)
    message("PCA variance plot saved to: ", save_plot)
  }

  # Print to console
  print(p)
  message("Top 10 PCs variance explained:")
  print(head(pca_df, 10))

  # Save Seurat object if path provided
  if (!is.null(save_path_rds)) {
    dir.create(dirname(save_path_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, file = save_path_rds)
    message("Seurat object with PCA saved to: ", save_path_rds)
  }

  # Store summary in misc slot
  seurat_obj@misc$pca_summary <- list(
    variance_df = pca_df,
    total_pcs = n_pcs
  )

  return(seurat_obj)
}
