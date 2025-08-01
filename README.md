# scRNAseqTools - Project Programming approaches for bioinformatics July 2025

This project provides a fully reproducible R package and Docker-based environment designed for the analysis of single-cell RNA sequencing data.

It guides users through the complete workflow - from gene annotation and quality filtering to dimensionality reduction, clustering, cell type identification, and tissue origin prediction - using Seurat and Bioconductor packages

---

## Docker Setup

The project includes a custom `Dockerfile` based on `rocker/rstudio:4.3.0` to ensure reproducibility.

## Dockerfile Summary

#### Dockerfile for R package development
```r
FROM rocker/rstudio:4.3.0
```

#### Install system dependencies required by R packages
```r
RUN apt-get update && apt-get install -y \
    libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev libgit2-dev \
    libudunits2-dev libgdal-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev && \
    rm -rf /var/lib/apt/lists/*
```

#### Install core CRAN packages for R package development
```r
RUN R -e "install.packages(c('devtools', 'usethis', 'roxygen2', 'testthat', 'knitr', 'rmarkdown', 'BiocManager'))"
```

#### Install Bioconductor packages for single-cell analysis
```r
RUN R -e "BiocManager::install(c('Seurat', 'SingleCellExperiment', 'scater', 'scran', 'SingleR', 'celldex', 'rtracklayer', 'BiocStyle'))"
```

#### Install supporting CRAN packages
```r
RUN R -e "install.packages(c('dplyr', 'ggplot2', 'plotly', 'viridis', 'cowplot', 'patchwork', 'Matrix', 'DT'))"
```

#### Expose RStudio Server on port 8787
```r
WORKDIR /home/rstudio/workspace
EXPOSE 8787
```

#### Build and run
```r
docker build -t scrnaseqtools .
docker run -e PASSWORD=yourpassword -p 8787:8787 -v ${PWD}:/home/rstudio/workspace scrnaseqtools
```

### Package structure
```r
scRNAseqTools/
├── R/                 # Modular R functions
├── man/               # Function documentation (.Rd)
├── vignettes/         # Full analysis workflow
├── DESCRIPTION, NAMESPACE
```
---

## Analysis Steps & Functions
Each analysis step is linked to a documented R function included in the package:
### Gene Annotation — filter_protein_coding_genes()
Filters raw gene counts to retain only protein-coding genes using a GTF annotation file.
```r
filter_protein_coding_genes <- function(seurat_obj, gtf_path) {
  gtf <- rtracklayer::import(gtf_path)
  gtf_genes <- gtf[gtf$type == "gene"]
  protein_coding_genes <- gtf_genes[gtf_genes$gene_type == "protein_coding"]
  coding_ids <- unique(protein_coding_genes$gene_id)
  seurat_obj <- subset(seurat_obj, features = intersect(rownames(seurat_obj), coding_ids))
  return(seurat_obj)
}
```
### Gene Expression Summary — plot_gene_expression_distribution()
Counts the number of genes with ≥3 UMIs per cell and displays a violin plot.
```r
plot_gene_expression_distribution <- function(seurat_obj) {
  umi_counts <- Matrix::colSums(seurat_obj@assays$RNA@counts >= 3)
  seurat_obj$UMI_3plus <- umi_counts
  VlnPlot(seurat_obj, features = "UMI_3plus") + 
    ggtitle("Genes with ≥3 UMIs per cell")
}
```
### Gene Filtering — filter_unwanted_genes()
Removes mitochondrial genes, ribosomal proteins, and ribosomal pseudogenes.
```r
filter_unwanted_genes <- function(seurat_obj) {
  mito_genes <- grep("^MT-", rownames(seurat_obj), value = TRUE)
  ribo_genes <- grep("^RP[SL]", rownames(seurat_obj), value = TRUE)
  pseudo_genes <- grep("RP[SL].*P[0-9]+$", rownames(seurat_obj), value = TRUE)

  all_remove <- unique(c(mito_genes, ribo_genes, pseudo_genes))
  summary_table <- data.frame(
    Category = c("Mitochondrial", "Ribosomal", "Ribosomal pseudogenes"),
    Removed = c(length(mito_genes), length(ribo_genes), length(pseudo_genes))
  )
  seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), all_remove))
  print(summary_table)
  return(seurat_obj)
}
```
#####                                    Category Genes_Removed
###### ribosomal_proteins       ribosomal_proteins            99
###### ribosomal_pseudogenes ribosomal_pseudogenes             0
###### mitochondrial_genes     mitochondrial_genes            13

### PCA — run_pca_analysis()
Performs PCA and plots the variance explained by the top 20 PCs.
```r
run_pca_analysis <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, npcs = 20)
    # Plot histogram
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC, y = Variance_Explained)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    ggplot2::labs(title = "Variance Explained by Principal Components",
                  x = "Principal Component",
                  y = "Proportion of Variance Explained") +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = 1:n_pcs)
  return(seuratobj)
}
```
### UMAP — run_umap_embedding()
Performs UMAP on selected PCs (e.g., first 15), and visualizes embedding.
```r
run_umap_embedding <- function(seurat_obj, dims = 1:15) {
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  DimPlot(seurat_obj, reduction = "umap", label = TRUE)
}
```
### Clustering — run_graph_clustering()
Runs graph-based clustering using Seurat.
```r
run_graph_clustering <- function(seurat_obj, dims = 1:15, resolution = 0.5) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
}
```
### Cell Type Annotation — annotate_cell_types()
Annotates cells using SingleR and HumanPrimaryCellAtlas reference.
```r
annotate_cell_types <- function(seurat_obj) {
  ref <- celldex::HumanPrimaryCellAtlasData()
  sce <- as.SingleCellExperiment(seurat_obj)
  pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
  seurat_obj$SingleR <- pred$labels
  DimPlot(seurat_obj, group.by = "SingleR", reduction = "umap", label = TRUE)
}
```
### Tissue Origin Inference
Interpretative step using known markers, cluster distributions, and reference labels. 
```r
markers <- Seurat::FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
 top_markers <- dplyr::group_by(markers, cluster) |>
    dplyr::top_n(n = 5, wt = avg_log2FC)
)
cell_type_summary <- table(seurat_obj$cell_type)
  print(cell_type_summary)
 dominant <- names(sort(cell_type_summary, decreasing = TRUE))[1:3]

  if (any(grepl("T_cell|B_cell|NK_cell", dominant))) {
    tissue_hypothesis <- "Lymphoid tissue (e.g., lymph node, spleen, thymus)"
  } else if (any(grepl("Monocyte|Macrophage|Dendritic", dominant))) {
    tissue_hypothesis <- "Immune-rich tissue (e.g., blood, bone marrow)"
  } else if (any(grepl("Fibroblast|Endothelial", dominant))) {
    tissue_hypothesis <- "Connective/vascular tissue"
  }
  message("Tissue origin hypothesis: ", tissue_hypothesis)
return(list(
    top_markers = top_markers,
    cell_type_summary = cell_type_summary,
    tissue_hypothesis = tissue_hypothesis
  ))
```
Based on marker gene expression, cluster structure and annotations. 

######           B_cell 853; CMP 4; GMP 4; HSC_-G-CSF 34; Macrophage 1; Monocyte 3807; Neutrophils 7; NK_cell 474; Platelets 27; Pre-Bcell_CD34 50; T_cells 4933
##### The inferred tissue of origin is: "Lymphoid tissue (e.g., lymph node, spleen, thymus)"

## Vignette
A complete documented pipeline, including all code, explanation, and results step by step, is provided in:
```r
vignettes/output/scrna_analysis_workflow.Rmd
```
And rendered as:
```r
vignettes/output/scrna_analysis_workflow.html
```

## Installation & Use
To install and use the package:
```r
devtools::document()
devtools::install()
library(scRNAseqTools)
```

## Docker Image
The complete analysis environment is available as a Docker image:
[gentilegiorgia01/scrnaseqtools](https://hub.docker.com/r/gentilegiorgia01/scrnaseqtools)

To run it locally:
```bash
docker pull gentilegiorgia01/scrnaseqtools:latest
docker run -e PASSWORD=yourpassword -p 8787:8787 gentilegiorgia01/scrnaseqtools
```
Then access RStudio in your browser at:
http://localhost:8787
Login -> 
Username: rstudio
Password: yourpassword
