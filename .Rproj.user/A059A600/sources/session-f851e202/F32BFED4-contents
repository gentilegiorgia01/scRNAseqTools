#' Filter Unwanted Genes from Seurat Object
#'
#' Removes mitochondrial genes, ribosomal proteins, and ribosomal pseudogenes
#' from a Seurat object based on a GTF annotation file.
#'
#' @param seurat_obj A Seurat object
#' @param gtf_path Path to the GTF file (e.g. Ensembl annotations)
#' @param save_path Optional path to save the filtered Seurat object (default: NULL)
#'
#' @return A filtered Seurat object with unwanted genes removed
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- filter_unwanted_genes(
#'   seurat_obj,
#'   gtf_path = "data/Homo_sapiens.GRCh38.111.gtf",
#'   save_path = "data/seurat_obj.filtered.rds"
#' )
#' }
#'
#' @importFrom rtracklayer import
filter_unwanted_genes <- function(seurat_obj, gtf_path, save_path = NULL) {

  # Validazioni iniziali
  if (!methods::is(seurat_obj, "Seurat")) stop("Input must be a Seurat object")
  if (!file.exists(gtf_path)) stop("GTF file not found at: ", gtf_path)

  message("Importing GTF file: ", gtf_path)
  gtf_data <- rtracklayer::import(gtf_path)
  gtf_genes <- gtf_data[gtf_data$type == "gene"]

  message("Filtering unwanted genes...")

  # Inizializza contatori
  genes_removed <- list(
    ribosomal_proteins = 0,
    ribosomal_pseudogenes = 0,
    mitochondrial_genes = 0
  )

  current_genes <- rownames(seurat_obj)

  # Ribosomi (RPS, RPL)
  ribo_genes <- grep("^RP[SL]", current_genes, value = TRUE)
  genes_removed$ribosomal_proteins <- length(ribo_genes)

  # Pseudogeni ribosomiali
  gene_match <- match(current_genes, gtf_genes$gene_name)
  biotypes <- gtf_genes$gene_biotype[gene_match]
  ribo_pseudo_genes <- current_genes[
    grepl("ribosomal", current_genes, ignore.case = TRUE) &
      grepl("pseudogene", biotypes, ignore.case = TRUE)
  ]
  genes_removed$ribosomal_pseudogenes <- length(ribo_pseudo_genes)

  # Geni mitocondriali (MT-)
  mito_genes <- grep("^MT-", current_genes, value = TRUE)
  genes_removed$mitochondrial_genes <- length(mito_genes)

  # Rimozione e aggiornamento Seurat object
  genes_to_remove <- unique(c(ribo_genes, ribo_pseudo_genes, mito_genes))
  genes_to_keep <- setdiff(current_genes, genes_to_remove)
  seurat_obj <- seurat_obj[genes_to_keep, ]

  # Riassunto rimozioni
  removal_summary <- data.frame(
    Category = names(genes_removed),
    Genes_Removed = unlist(genes_removed),
    stringsAsFactors = FALSE
  )

  message("Gene removal summary:")
  print(removal_summary)

  message("Remaining genes after filtering: ", nrow(seurat_obj))

  # Salvataggio opzionale
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(seurat_obj, save_path)
    message("Filtered object saved to: ", save_path)
  }

  # Salva riassunto nei metadati
  seurat_obj@misc$gene_filtering <- list(
    summary = removal_summary,
    gtf_file = gtf_path,
    genes_remaining = nrow(seurat_obj)
  )

  return(seurat_obj)
}
