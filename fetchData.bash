#!/bin/bash

# Robjects
echo "Fetching R objects..."
synapse get syn21560469 --downloadLocation ./Analysis/seurat/
synapse get syn21560492 --downloadLocation ./Analysis/seurat/
synapse get syn21560508 --downloadLocation ./Analysis/seurat/
synapse get syn21560419 --downloadLocation ./Analysis/seurat/
synapse get syn21560426 --downloadLocation ./Analysis/seurat/
synapse get syn21624617 --downloadLocation ./Analysis/seurat/
synapse get syn21624978 --downloadLocation ./Analysis/seurat/
synapse get syn21560428 --downloadLocation ./Analysis/seurat/

echo "Fetching cellphonedb tsv..."
synapse get syn21560548 --downloadLocation ./Analysis/cellphonedb/
