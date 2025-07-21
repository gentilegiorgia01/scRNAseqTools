#' Annotate Protein-Coding Genes
#'
#' Identifies and filters protein-coding genes from GTF annotation file and
#' retains only these genes in the Seurat object.
#'
#' @param data_path Character. Path to the single-cell dataset (H5 format)
#' @param gtf_path Character. Path to the GTF annotation file
#' @param project_name Character. Name for the Seurat project (default: "scRNA_analysis")
#'
#' @return A Seurat object containing only protein-coding genes
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_obj <- annotate_protein_coding(
#'   data_path = "data/filtered_feature_bc_matrix",
#'   gtf_path = "data/Homo_sapiens.GRCh38.111.gtf"
#' )
#' }
#'
#' @importFrom Seurat Read10X_h5 CreateSeuratObject
#' @importFrom rtracklayer import
#' @importFrom methods is
annotate_protein_coding <- function(data_path, gtf_path, project_name = "scRNA_analysis") {

  # Input validation
  if (!file.exists(data_path)) {
    stop("Data file not found: ", data_path)
  }
  if (!file.exists(gtf_path)) {
    stop("GTF file not found: ", gtf_path)
  }

  # Load GTF file and extract protein-coding genes
  message("Loading GTF file...")
  gtf_data <- rtracklayer::import(gtf_path)

  # Filter for protein-coding genes
  protein_coding_genes <- gtf_data[gtf_data$gene_biotype == "protein_coding"]
  protein_coding_gene_names <- unique(protein_coding_genes$gene_name)

  message("Found ", length(protein_coding_gene_names), " protein-coding genes")

  # Load single-cell data
  message("Loading single-cell data...")
  count_data <- Seurat::Read10X(data_path)
  seurat_obj <- Seurat::CreateSeuratObject(counts = count_data, project = project_name)

  # Filter for protein-coding genes only
  genes_to_keep <- intersect(rownames(seurat_obj), protein_coding_gene_names)
  seurat_obj <- seurat_obj[genes_to_keep, ]

  message("Retained ", length(genes_to_keep), " protein-coding genes in the dataset")

  # Store annotation info in misc slot
  seurat_obj@misc$annotation_info <- list(
    total_protein_coding_genes = length(protein_coding_gene_names),
    genes_in_dataset = length(genes_to_keep),
    gtf_file = gtf_path
  )

  return(seurat_obj)
}
