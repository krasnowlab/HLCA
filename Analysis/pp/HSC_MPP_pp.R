# Read in MPP and HSC expression profiles
MPP <- read.csv(here('databulk', 'MPP1_CGTACTAG_TATCCTCT.2p.ReadsPerGene.out.tab'), sep = "\t", row.names = 1)
HSC <- read.csv(here('databulk', 'HSC1_TAAGGCGA_CTCTCTAT.2p.ReadsPerGene.out.tab'), sep = "\t", row.names = 1)
bulkSeq2 <- cbind(MPP[,1],HSC[,1])
rownames(bulkSeq2) <- rownames(MPP)
bulkSeq2 <- bulkSeq2[-grep('^N_', rownames(bulkSeq2)),]
bulkSeq2 <- as.data.frame(bulkSeq2)

# convert genome ids to hgnc symbols
humanGTF <- read.csv('../Demux/gencode.v29.chr_patch_hapl_scaff.annotation_w_ERCC.gtf', sep = "\t", quote = "", header = FALSE)
humanGTF <- humanGTF[,'V9']
humanGTF <- humanGTF[1:2989176]

geneIDs <- gsub('.*?gene_id "([^"]+)";.*', '\\1', humanGTF)
geneSyms <- gsub('.*?gene_name "([^"]+)";.*', '\\1', humanGTF)
geneSyms[which(geneSyms == humanGTF)] = geneIDs[which(geneSyms == humanGTF)]

geneMatrix <- matrix(0,nrow = 2989176, ncol = 2)
geneMatrix[,1] <- geneIDs
geneMatrix[,2] <- geneSyms
geneMatrix <- geneMatrix[-intersect(which(duplicated(geneIDs)), which(duplicated(geneSyms))),]
geneMatrix <- as.data.frame(geneMatrix)
colnames(geneMatrix) <- c('geneID','geneSym')
rownames(geneMatrix) <- geneMatrix[,'geneID']

toSet <- geneMatrix[rownames(bulkSeq2),'geneSym']
toSet <- as.character(toSet)
toSet[which(is.na(toSet))] <- rownames(bulkSeq2)[which(is.na(toSet))]
bulkSeq2[,'geneSym'] <- toSet


one2one <- bulkSeq2[-which(toSet %in% toSet[duplicated(toSet)]),]
rownames(one2one) <- one2one[,'geneSym']
one2one[,'geneSym'] <- NULL

complex <- bulkSeq2[which(toSet %in% toSet[duplicated(toSet)]),]
complex <- aggregate(. ~ geneSym, complex, sum)
rownames(complex) <- complex[,'geneSym']
complex[,'geneSym'] <- NULL

bulkSeq2 <- rbind(one2one, complex)
colnames(bulkSeq2) <- c('MPP', 'HSC')
bulkSeq2 <- log1p(bulkSeq2)

bulkSeq2 <- bulkSeq2[intersect(rownames(bulkSeq2), rownames(bulkSeq)),]
bulkSeq <- bulkSeq[intersect(rownames(bulkSeq2), rownames(bulkSeq)),]
