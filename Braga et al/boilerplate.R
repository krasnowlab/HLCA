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

load_tissue_droplet = function(dataset){
  
  raw.data <- read.csv(file = here('datadropseq', paste(dataset, 'raw_counts.csv', sep = "_")), row.names = 1)
  meta.data <- read.csv(file = here('metadata', paste(dataset, 'barcodes_cell_types.txt', sep="_")), row.names = 1, sep = "\t")
  
  # Create the Seurat object with all the data
  tiss <- CreateSeuratObject(raw.data = raw.data, project = 'Barga et al')
  tiss <- AddMetaData(object = tiss, metadata = meta.data)
  
  # Create metadata columns for annotations
  tiss@meta.data[,'free_annotation'] <- NA
  
  tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nUMI"), 
                      low.thresholds = c(100, 500))
  tiss <- process_tissue(tiss, 1e4)
  
  return(tiss)
}