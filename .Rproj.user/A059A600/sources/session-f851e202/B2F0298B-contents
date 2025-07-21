#' Annotate Cell Types Using SingleR and HumanPrimaryCellAtlas
#'
#' Annotates cells in a Seurat object using SingleR and a reference dataset,
#' and plots the results on the UMAP embedding.
#'
#' @param seurat_obj A Seurat object with PCA and UMAP already computed
#' @param save_path_plot Character. Path to save the UMAP cell type plot (default: NULL)
#' @param save_path_csv Character. Path to save the cluster vs cell type table (default: NULL)
#' @param save_path_rds Character. Path to save the updated Seurat object (default: NULL)
#'
#' @return Seurat object with cell type annotations added
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- annotate_cell_types_singleR(
#'   seurat_obj,
#'   save_path_plot = "data/umap_celltype_annotations.png",
#'   save_path_csv = "data/cluster_vs_celltype.csv",
#'   save_path_rds = "data/seurat_obj.annotated.rds"
#' )
#' }
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SingleR SingleR
#' @importFrom celldex HumanPrimaryCellAtlasData
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 ggsave ggtitle theme_minimal theme
annotate_cell_types_singleR <- function(seurat_obj,
                                        save_path_plot = NULL,
                                        save_path_csv = NULL,
                                        save_path_rds = NULL) {
  if (!methods::is(seurat_obj, "Seurat")) stop("Input must be a Seurat object")

  message("Loading reference dataset...")
  ref_data <- celldex::HumanPrimaryCellAtlasData()

  message("Converting to SingleCellExperiment...")
  sce <- as.SingleCellExperiment(seurat_obj)

  message("Running SingleR annotation...")
  annotations <- SingleR::SingleR(test = sce, ref = ref_data, labels = ref_data$label.main)
  seurat_obj$cell_type <- annotations$labels

  message("Plotting UMAP with cell type annotations...")
  annotation_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type") +
    ggplot2::ggtitle("Cell Type Annotations") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  print(annotation_plot)

  if (!is.null(save_path_plot)) {
    dir.create(dirname(save_path_plot), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(save_path_plot, plot = annotation_plot, width = 7, height = 5)
    message("Annotation plot saved to: ", save_path_plot)
  }

  if (!is.null(save_path_csv)) {
    comparison_table <- table(seurat_obj$seurat_clusters, seurat_obj$cell_type)
    write.csv(as.data.frame(comparison_table), save_path_csv, row.names = TRUE)
    message("Cluster vs Cell Type table saved to: ", save_path_csv)
  }

  if (!is.null(save_path_rds)) {
    saveRDS(seurat_obj, save_path_rds)
    message("Annotated Seurat object saved to: ", save_path_rds)
  }

  # Add metadata
  seurat_obj@misc$annotation_info <- list(
    reference = "HumanPrimaryCellAtlas",
    method = "SingleR",
    unique_cell_types = length(unique(seurat_obj$cell_type))
  )

  return(seurat_obj)
}
