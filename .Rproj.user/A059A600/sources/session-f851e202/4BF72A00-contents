#' Perform Clustering and Plot UMAP with Clusters
#'
#' Executes clustering on the Seurat object using the Louvain algorithm and visualizes the result.
#'
#' @param seurat_obj A Seurat object with UMAP and PCA already computed
#' @param n_pcs Integer. Number of PCs to use for clustering (default: 15)
#' @param resolution Numeric. Clustering resolution parameter (default: 0.5)
#' @param algorithm Integer. Clustering algorithm (1 = Louvain, 2 = Louvain multilevel, 3 = SLM) (default: 1)
#' @param save_path_plot Character. Path to save the cluster UMAP plot (default: NULL)
#' @param save_path_rds Character. Path to save the updated Seurat object (default: NULL)
#'
#' @return Seurat object with clustering information added
#' @export
run_clustering_and_plot <- function(seurat_obj,
                                    n_pcs = 15,
                                    resolution = 0.5,
                                    algorithm = 1,
                                    save_path_plot = NULL,
                                    save_path_rds = NULL) {

  if (!methods::is(seurat_obj, "Seurat")) stop("Input must be a Seurat object")

  message("Performing clustering using ", n_pcs, " PCs, resolution ", resolution, ", algorithm ", algorithm, "...")

  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:n_pcs)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution, algorithm = algorithm)

  cluster_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
    ggplot2::ggtitle("Clustering Results") +
    ggplot2::theme_minimal()

  print(cluster_plot)

  if (!is.null(save_path_plot)) {
    dir.create(dirname(save_path_plot), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path_plot, plot = cluster_plot, width = 6, height = 5)
    message("Clustering plot saved to: ", save_path_plot)
  }

  if (!is.null(save_path_rds)) {
    dir.create(dirname(save_path_rds), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, save_path_rds)
    message("Seurat object with clustering saved to: ", save_path_rds)
  }

  seurat_obj@misc$clustering_info <- list(
    n_pcs = n_pcs,
    resolution = resolution,
    algorithm = algorithm,
    n_clusters = length(unique(seurat_obj$seurat_clusters))
  )

  return(seurat_obj)
}
