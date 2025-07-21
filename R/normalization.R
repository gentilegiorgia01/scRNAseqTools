#' Perform optimized normalization for large datasets
#'
#' @param seurat_obj Seurat object
#' @param method Normalization method ("LogNormalize" or "SCT")
#' @param verbose Logical, show progress messages
#' @param gc_freq Frequency of garbage collection (every n steps)
#'
#' @return Normalized Seurat object
#' @export
perform_optimized_normalization <- function(seurat_obj,
                                            method = "LogNormalize",
                                            verbose = TRUE,
                                            gc_freq = 1000) {

  if (verbose) {
    cat("Starting optimized normalization...\n")
    cat("Method:", method, "\n")
    cat("Cells:", ncol(seurat_obj), "\n")
    cat("Genes:", nrow(seurat_obj), "\n")
    cat("Memory before normalization:\n")
    print(gc())
  }

  # Forza garbage collection prima di iniziare
  invisible(gc(verbose = FALSE))

  if (method == "LogNormalize") {
    # Normalizzazione standard ottimizzata
    if (verbose) cat("Performing LogNormalize...\n")

    seurat_obj <- Seurat::NormalizeData(
      seurat_obj,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = verbose
    )

  } else if (method == "SCT") {
    # SCTransform per dataset più grandi
    if (verbose) cat("Performing SCTransform...\n")

    seurat_obj <- Seurat::SCTransform(
      seurat_obj,
      verbose = verbose,
      return.only.var.genes = TRUE  # Riduce memoria
    )
  }

  # Garbage collection dopo normalizzazione
  invisible(gc(verbose = FALSE))

  if (verbose) {
    cat("Normalization completed!\n")
    cat("Memory after normalization:\n")
    print(gc())
  }

  return(seurat_obj)
}

#' Perform PCA with memory optimization
#'
#' @param seurat_obj Seurat object (normalized)
#' @param npcs Number of principal components to compute
#' @param verbose Logical, show progress
#'
#' @return Seurat object with PCA computed
#' @export
perform_optimized_pca <- function(seurat_obj, npcs = 50, verbose = TRUE) {

  if (verbose) {
    cat("Finding variable features...\n")
  }

  # Trova feature variabili se non già fatto
  if (!"var.features" %in% names(seurat_obj@assays[[Seurat::DefaultAssay(seurat_obj)]]@meta.features) ||
      length(Seurat::VariableFeatures(seurat_obj)) == 0) {

    seurat_obj <- Seurat::FindVariableFeatures(
      seurat_obj,
      selection.method = "vst",
      nfeatures = 2000,
      verbose = verbose
    )
  }

  if (verbose) {
    cat("Variable features found:", length(Seurat::VariableFeatures(seurat_obj)), "\n")
    cat("Scaling data...\n")
  }

  # Scala solo le feature variabili per risparmiare memoria
  seurat_obj <- Seurat::ScaleData(
    seurat_obj,
    features = Seurat::VariableFeatures(seurat_obj),
    verbose = verbose
  )

  # Garbage collection
  invisible(gc(verbose = FALSE))

  if (verbose) {
    cat("Running PCA...\n")
  }

  # PCA
  seurat_obj <- Seurat::RunPCA(
    seurat_obj,
    features = Seurat::VariableFeatures(seurat_obj),
    npcs = npcs,
    verbose = verbose
  )

  # Final garbage collection
  invisible(gc(verbose = FALSE))

  if (verbose) {
    cat("PCA completed!\n")
    print(gc())
  }

  return(seurat_obj)
}
