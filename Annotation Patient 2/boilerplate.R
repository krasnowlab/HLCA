library(useful)
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)
my.cols <- rev(brewer.pal(11, "RdYlBu"))

load_tissue_facs = function(tissue_of_interest){
  
  # Load the per-plate metadata
  plate_metadata_filename = here('metadata', 'SS2.csv')

  plate_metadata <- read.csv(plate_metadata_filename, sep=",", header = TRUE)
  colnames(plate_metadata)[1] <- "plate.barcode"
  
  
  # Load the gene names and set the metadata columns by opening the first file
  filename = here('datass2', paste0(tissue_of_interest, '-counts.csv'))
  
  raw.data = read.csv(filename, sep=",", row.names=1)
  
  plate.barcodes <- unlist(lapply(colnames(raw.data), function(x) strsplit(x, "\\.")[[1]][1]))
  plate.barcodes = lapply(plate.barcodes, function(x) strsplit(x, "_")[[1]][2])
  
  barcode.df = t.data.frame(as.data.frame(plate.barcodes))
  
  rownames(barcode.df) = colnames(raw.data)
  colnames(barcode.df) = c('plate.barcode')
  barcode.df = as.data.frame(barcode.df)

  barcode.df[,'cell.id'] = rownames(barcode.df)
  meta.data <- merge(barcode.df, plate_metadata, by='plate.barcode', sort = F)
  rownames(meta.data) <- meta.data[,'cell.id']
  
  # Sort cells by cell name
  meta.data = meta.data[order(rownames(meta.data)), ]
  raw.data = raw.data[, rownames(meta.data)]
  
  raw.data = raw.data[-grep("N_", rownames(raw.data)),]
  
  # Find ERCC's, compute the percent ERCC, and drop them from the raw data.
  erccs <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = TRUE)
  percent.ercc <- Matrix::colSums(raw.data[erccs, ])/Matrix::colSums(raw.data)
  ercc.index <- grep(pattern = "^ERCC-", x = rownames(x = raw.data), value = FALSE)
  ercc.matrix <- as.data.frame(t(raw.data[ercc.index,]))
  raw.data <- raw.data[-ercc.index,]
  
  # Load the STAR alignment stats
  filename = here('patient_stats', '2.csv')
  star.data = as.data.frame(t(read.csv(filename, sep=",", row.names=1)))
  
  # Create the Seurat object with all the data
  tiss <- CreateSeuratObject(raw.data = raw.data, project = "Human Lung SS2 - Patient 2")
  tiss <- AddMetaData(object = tiss, meta.data)
  tiss <- AddMetaData(object = tiss, percent.ercc, col.name = "percent.ercc")
  
  
  # Change default name for sums of counts from nUMI to nReads
  colnames(tiss@meta.data)[colnames(tiss@meta.data) == 'nUMI'] <- 'nReads'
  
  
  #Calculate percent ribosomal genes.
  ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = tiss@data), value = TRUE)
  percent.ribo <- Matrix::colSums(tiss@raw.data[ribo.genes, ])/Matrix::colSums(tiss@raw.data)
  tiss <- AddMetaData(object = tiss, metadata = percent.ribo, col.name = "percent.ribo")
  
  
  # Create metadata columns for cell annotations
  tiss@meta.data[,'gating'] <- NA
  tiss@meta.data[,'free_annotation'] <- NA
  
  tiss <- AddMetaData(object = tiss, star.data)
  tiss <- AddMetaData(object = tiss, ercc.matrix)
  
  tmp <- colnames(tiss@meta.data)
  
  tiss <- FilterCells(object = tiss, subset.names = c("nGene", "nReads"), 
                      low.thresholds = c(500, 50000))
  
  colnames(tiss@meta.data) <- tmp
  
  tiss <- process_tissue(tiss, 1e6)
  
  return(tiss)
}

process_tissue = function(tiss, scale){
  tiss <- NormalizeData(object = tiss, scale.factor = scale)
  tiss <- ScaleData(object = tiss)
  tiss <- FindVariableGenes(object = tiss, do.plot = TRUE, x.high.cutoff = Inf, y.cutoff = 0.5)
  tiss <- RunPCA(object = tiss, do.print = FALSE)
  tiss <- ProjectPCA(object = tiss, do.print = FALSE)
}

load_tissue_droplet = function(tissue_of_interest){
  
  # read the metadata to get the channels we want
  droplet_metadata_filename = here('metadata', '10X.csv')
  
  droplet_metadata <- read.csv(droplet_metadata_filename, sep=",", header = TRUE)
  colnames(droplet_metadata)[1] <- "channel"
  
  if (tissue_of_interest != 'all') {
    tissue_metadata = filter(droplet_metadata, tissue == tissue_of_interest)[,c('channel','tissue','region','compartment','sample','location','patient')]
  } else {
    tissue_metadata = droplet_metadata[,c('channel','tissue','region','compartment','sample','location','patient')]
  }
  
  subfolder = tissue_metadata$channel[1]
  raw.data <- Read10X(data.dir = here('data10x', subfolder))
  colnames(raw.data) <- lapply(colnames(raw.data), function(x) paste0(tissue_metadata$channel[1], '_', x))
  meta.data = data.frame(row.names = colnames(raw.data))
  meta.data['channel'] = tissue_metadata$channel[1]
  
  if (length(tissue_metadata$channel) > 1){
    # Some tissues, like Thymus and Heart had only one channel
    for(i in 2:nrow(tissue_metadata)){
      subfolder = tissue_metadata$channel[i]
      new.data <- Read10X(data.dir = here('data10x', subfolder))
      colnames(new.data) <- lapply(colnames(new.data), function(x) paste0(tissue_metadata$channel[i], '_', x))
      
      new.metadata = data.frame(row.names = colnames(new.data))
      new.metadata['channel'] = tissue_metadata$channel[i]
      
      raw.data = cbind(raw.data, new.data)
      meta.data = rbind(meta.data, new.metadata)
    }
  }
  
  rnames = row.names(meta.data)
  meta.data <- merge(meta.data, tissue_metadata, sort = F)
  row.names(meta.data) <- rnames
  
  # Order the cells alphabetically to ensure consistency.
  
  ordered_cell_names = order(colnames(raw.data))
  raw.data = raw.data[,ordered_cell_names]
  meta.data = meta.data[ordered_cell_names,]
  
  # Create the Seurat object with all the data
  tiss <- CreateSeuratObject(raw.data = raw.data, project = "Human Lung 10x - Patient 2")
  
  tiss <- AddMetaData(object = tiss, meta.data)

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