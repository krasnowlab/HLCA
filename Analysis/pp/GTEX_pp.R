GTEX_metadata <- read.csv(file = here('databulk', 'GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt'), row.names = 1, sep = "\t")

HL_GTEX_metadata <- GTEX_metadata[setdiff(intersect(intersect(intersect(which(GTEX_metadata$SMAFRZE %in% c('RNASEQ', 'OMNI')),which(GTEX_metadata$SMTS == 'Lung')),which(GTEX_metadata$SMATSSCR < 2)),which(GTEX_metadata$SMTSISCH < 250)), union(grep('(moderate|severe|dense|chronic|marked|some|damage|patchy|acute|dilated|failure)', GTEX_metadata$SMPTHNTS, ignore.case = TRUE), grep('(emphys|pneumonia|granulomas|atelectasis|pneumonitis|atelectatic)', GTEX_metadata$SMPTHNTS, ignore.case = TRUE))),]

rownames(HL_GTEX_metadata) <- gsub('-', '.', rownames(HL_GTEX_metadata))

#

GTEX_data <- read.csv(file = here('databulk', 'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.txt'), row.names = 1, sep = "\t")
HL_GTEX_data <- GTEX_data[,c('Description', rownames(HL_GTEX_metadata))]

HL_GTEX_data <- HL_GTEX_data %>% group_by(Description) %>% summarise_all(funs(sum))

HL_GTEX_data <- as.data.frame(HL_GTEX_data)

rownames(HL_GTEX_data) <- HL_GTEX_data$Description

HL_GTEX_data[,'Description'] <- NULL
