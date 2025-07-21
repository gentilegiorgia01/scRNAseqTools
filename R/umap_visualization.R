#' Run UMAP and Generate UMAP Plot
#'
#' Runs UMAP on the given Seurat object using the specified number of PCs and plots the result.
#'
#' @param seurat_obj A Seurat object with PCA already computed
#' @param n_pcs Integer. Number of principal components to use for UMAP (default: 15)
#' @param save_path_plot Character. Path to save the UMAP plot (default: NULL)
#' @param save_path_rds Character. Path to save the updated Seurat object (default: NULL)
#'
#' @return A Seurat object with UMAP embedding added
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- run_umap_and_plot(
#'   seurat_obj,
#'   n_pcs = 15,
#'   save_path_plot = "data/umap_plot.png",
#'   save_path_rds = "data/seurat_obj.umap.rds"
#' )
#' }
#'
#' @importFrom Seurat RunUMAP DimPlot
#' @importFrom ggplot2 ggsave ggtitle theme_minimal
run_umap_and_plot <- function(seurat_obj,
                              n_pcs = 15,
                              save_path_plot = NULL,
                              save_path_rds = NULL) {

  if (!methods::is(seurat_obj, "Seurat")) stop("Input must be a Seurat object")
  if (!is.numeric(n_pcs) || n_pcs <= 0) stop("n_pcs must be a positive integer")

  message("Running UMAP using the first ", n_pcs, " principal components...")

  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:n_pcs)

  umap_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap") +
    ggplot2::ggtitle("UMAP Visualization") +
    ggplot2::theme_minimal()

  print(umap_plot)

  if (!is.null(save_path_plot)) {
    dir.create(dirname(save_path_plot), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path_plot, plot = umap_plot, width = 6, height = 5)
    message("UMAP plot saved to: ", save_path_plot)
  }

  if (!is.null(save_path_rds)) {
    dir.create(dirname(save_path_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, save_path_rds)
    message("Seurat object with UMAP saved to: ", save_path_rds)
  }

  seurat_obj@misc$umap_info <- list(
    dims_used = n_pcs
  )

  return(seurat_obj)
}
