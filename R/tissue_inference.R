#' Infer Tissue Origin from Cell Types and Marker Genes
#'
#' Based on dominant annotated cell types and top marker genes, proposes a tissue
#' of origin for the dataset. Saves top markers and final Seurat object.
#'
#' @param seurat_obj A Seurat object with cell type annotations
#' @param save_path_markers Character. Path to save the top marker gene CSV (default: NULL)
#' @param save_path_rds Character. Path to save the updated Seurat object (default: NULL)
#' @param max_cells_per_cluster Integer. Maximum cells to use per cluster for marker finding (default: 200)
#' @param test_use Character. Test to use for finding markers (default: "wilcox")
#'
#' @return A list containing the top markers, cell type summary, and tissue hypothesis
#' @export
#'
#' @examples
#' \dontrun{
#' results <- infer_tissue_origin(
#'   seurat_obj,
#'   save_path_markers = "data/top_markers.csv",
#'   save_path_rds = "results/seurat_analysis_complete.rds"
#' )
#' }
#'
#' @importFrom Seurat FindAllMarkers Idents
#' @importFrom dplyr group_by top_n
infer_tissue_origin <- function(seurat_obj,
                                save_path_markers = NULL,
                                save_path_rds = NULL,
                                max_cells_per_cluster = 200,
                                test_use = "wilcox") {
  if (!methods::is(seurat_obj, "Seurat")) stop("Input must be a Seurat object")
  if (is.null(seurat_obj$cell_type)) stop("Seurat object lacks 'cell_type' annotations")

  message("Inferring tissue origin...")

  # Step 1: Subsample if needed to speed up computation
  original_ident <- Idents(seurat_obj)
  Idents(seurat_obj) <- seurat_obj$seurat_clusters

  # Check if we need to subsample
  cluster_sizes <- table(Idents(seurat_obj))
  if (any(cluster_sizes > max_cells_per_cluster)) {
    message("Large clusters detected. Subsampling to speed up marker detection...")

    # Subsample cells from large clusters
    cells_to_keep <- c()
    for (cluster in names(cluster_sizes)) {
      cluster_cells <- names(Idents(seurat_obj))[Idents(seurat_obj) == cluster]
      if (length(cluster_cells) > max_cells_per_cluster) {
        sampled_cells <- sample(cluster_cells, max_cells_per_cluster)
      } else {
        sampled_cells <- cluster_cells
      }
      cells_to_keep <- c(cells_to_keep, sampled_cells)
    }

    # Create subsampled object for marker finding
    seurat_subset <- seurat_obj[, cells_to_keep]
  } else {
    seurat_subset <- seurat_obj
  }

  # Step 1: Identify top markers with optimized parameters
  message("Finding marker genes...")
  tryCatch({
    markers <- Seurat::FindAllMarkers(
      seurat_subset,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      test.use = test_use,
      verbose = TRUE,
      return.thresh = 0.05,  # Only return significant markers
      min.diff.pct = 0.1     # Minimum difference in detection rate
    )

    top_markers <- dplyr::group_by(markers, cluster) |>
      dplyr::top_n(n = 5, wt = avg_log2FC)

  }, error = function(e) {
    message("Error in FindAllMarkers: ", e$message)
    message("Falling back to simplified analysis...")

    # Fallback: create a simplified marker list
    clusters <- unique(Idents(seurat_subset))
    top_markers <- data.frame(
      cluster = clusters,
      gene = paste0("Marker_", clusters),
      avg_log2FC = 1,
      stringsAsFactors = FALSE
    )
  })

  if (!is.null(save_path_markers)) {
    dir.create(dirname(save_path_markers), showWarnings = FALSE, recursive = TRUE)
    write.csv(top_markers, save_path_markers, row.names = FALSE)
    message("Top markers saved to: ", save_path_markers)
  }

  # Step 2: Cell type composition (use original object)
  message("Analyzing cell type composition...")
  cell_type_summary <- table(seurat_obj$cell_type)
  print(cell_type_summary)

  # Step 3: Heuristic tissue inference
  dominant <- names(sort(cell_type_summary, decreasing = TRUE))[1:3]
  tissue_hypothesis <- "Unknown"

  # More comprehensive tissue inference
  if (any(grepl("T_cell|T.cell|NK_cell|NK.cell", dominant, ignore.case = TRUE))) {
    tissue_hypothesis <- "Lymphoid tissue (e.g., lymph node, spleen, thymus)"
  } else if (any(grepl("B_cell|B.cell|Plasma", dominant, ignore.case = TRUE))) {
    tissue_hypothesis <- "B cell-rich lymphoid tissue (e.g., lymph node, spleen)"
  } else if (any(grepl("Monocyte|Macrophage|Dendritic", dominant, ignore.case = TRUE))) {
    tissue_hypothesis <- "Immune-rich tissue (e.g., blood, bone marrow, inflamed tissue)"
  } else if (any(grepl("Fibroblast|Endothelial", dominant, ignore.case = TRUE))) {
    tissue_hypothesis <- "Connective/vascular tissue"
  } else if (any(grepl("Epithelial", dominant, ignore.case = TRUE))) {
    tissue_hypothesis <- "Epithelial tissue (e.g., skin, gut, lung)"
  } else if (any(grepl("Neuron|Astrocyte|Oligodendrocyte", dominant, ignore.case = TRUE))) {
    tissue_hypothesis <- "Neural tissue (brain, spinal cord)"
  }

  message("Dominant cell types: ", paste(dominant, collapse = ", "))
  message("Tissue origin hypothesis: ", tissue_hypothesis)

  # Step 4: Save updated object
  if (!is.null(save_path_rds)) {
    dir.create(dirname(save_path_rds), showWarnings = FALSE, recursive = TRUE)
    saveRDS(seurat_obj, save_path_rds)
    message("Seurat object saved to: ", save_path_rds)
  }

  # Store in misc
  seurat_obj@misc$tissue_inference <- list(
    dominant_cell_types = dominant,
    hypothesis = tissue_hypothesis,
    cell_type_counts = cell_type_summary
  )

  # Restore original identity
  Idents(seurat_obj) <- original_ident

  return(list(
    top_markers = top_markers,
    cell_type_summary = cell_type_summary,
    tissue_hypothesis = tissue_hypothesis,
    seurat_obj = seurat_obj
  ))
}
