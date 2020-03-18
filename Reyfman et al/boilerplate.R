if (!require("useful")) {
  install.packages("useful", dependencies = TRUE)
  library(useful)
}
if (!require("here")) {
  install.packages("here", dependencies = TRUE)
  library(here)
}
if (!require("Seurat")) {
  install.packages("Seurat", dependencies = TRUE)
  library(Seurat)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("Matrix")) {
  install.packages("Matrix", dependencies = TRUE)
  library(Matrix)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
my.cols <- rev(brewer.pal(11, "RdYlBu"))

process_tissue = function(tiss, scale){
  tiss <- NormalizeData(object = tiss, scale.factor = scale)
  tiss <- ScaleData(object = tiss)
  tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)
  tiss <- RunPCA(object = tiss, do.print = FALSE)
  tiss <- ProjectPCA(object = tiss, do.print = FALSE)
}

load_tissue_droplet = function(datasets){
  
  for (dataset in datasets) {
    new.data <- Read10X_h5(filename = here('data10x', paste0(dataset, '_filtered_gene_bc_matrices_h5.h5')))
    colnames(new.data) <- paste(dataset, colnames(new.data), sep = "_")
    if (exists('raw.data')) {
      raw.data <- cbind(raw.data, new.data)
    } else {
      raw.data <- new.data
    } 
  }
  
  # Create the Seurat object with all the data
  tiss <- CreateSeuratObject(raw.data = raw.data, project = 'Reyfman et al')
  
  ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = tiss@data), value = TRUE)
  percent.ribo <- Matrix::colSums(tiss@raw.data[ribo.genes, ])/Matrix::colSums(tiss@raw.data)
  tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")
  
  # Create metadata columns for annotations
  tiss@meta.data[,'free_annotation'] <- NA
  
  tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nUMI"), 
                      low.thresholds = c(500, 1000))
  tiss <- process_tissue(tiss, 1e4)
  
  return(tiss)
}