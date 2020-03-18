PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}
CDotPlot <- function (object, genes.plot, cols.use = c("lightgrey", "blue"), 
          col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, scale.min = NA, scale.max = NA,
          group.by, plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE, levels.order = FALSE, only.exp = FALSE, ident.include = levels(object@ident), dot.outline = 'white', flip.axis = FALSE, pct.range = c(0.1,1)) 
{
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  object@ident <- factor(object@ident, levels = levels(object@ident)[levels.order])
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot, cells.use = WhichCells(object, ident = ident.include)))
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident[WhichCells(object, ident = ident.include)]
  data.to.plot <- data.to.plot %>% gather(key = genes.plot, value = expression, -c(cell, id))
  if (only.exp == TRUE) {
    data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
      summarize(avg.exp = mean(log2(expm1(x = expression[expression > 0]) + 1)), pct.exp = PercentAbove(x = expression, threshold = 0))
  } else {
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
    summarize(avg.exp = mean(log2(expm1(x = expression) + 1)), pct.exp = PercentAbove(x = expression, threshold = 0))
  }
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    mutate(avg.exp = MinMax(data = avg.exp, max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                    levels = sub(pattern = "-", replacement = ".", 
                                                         x = genes.plot), ordered = TRUE)

  data.to.plot$pct.exp[data.to.plot$pct.exp <= dot.min] <- NA
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                 y = id)) + geom_point(mapping = aes(size = pct.exp, 
                                                                                     fill = avg.exp), colour=dot.outline,pch=21) + scale_radius(range = c(0,
                                                                                                                                                      dot.scale), limits = pct.range) + scale_fill_gradientn(colors = cols.use, guide = guide_colorbar(ticks = FALSE)) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, 
                                            hjust = 1))
  }
  if (!flip.axis) {
    p <- p + coord_flip()
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}

NEDotPlot <- function(barcodes, cols.use = c("lightgrey", "blue"), 
          col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, scale.min = NA, scale.max = NA,
          group.by, plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE, levels.order = FALSE, only.exp = FALSE) 
{
  barcodes <- as.matrix(barcodes)
  barcodes <- barcodes > 0
  barcodes[which(barcodes == TRUE)] <- 1
  barcodes <- as.data.frame(t(barcodes))
  data.to.plot <- barcodes
  ids <- c()
  counts <- numeric()
  for (i in 1:dim(barcodes)[1]) {
    ids <- c(ids,paste(as.numeric(barcodes[i,]), collapse=""))
    counts <- c(counts,sum(as.numeric(barcodes[i,])))
    
  }

  names(ids) <- rownames(barcodes)
  ids <- as.factor(ids)
  if (levels.order == FALSE) {
    ids <- factor(ids, levels = levels(ids)[order(sapply(levels(ids), function(x) sum(as.numeric(strsplit(as.character(x), split='')[[1]]))), decreasing = TRUE)])
  }
  genes.plot <- colnames(barcodes)
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- ids
  data.to.plot <- data.to.plot %>% gather(key = genes.plot, value = expression, -c(cell, id))
  if (only.exp == TRUE) {
    data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
      summarize(avg.exp = mean(log2(expm1(x = expression[expression > 0]) + 1)), pct.exp = PercentAbove(x = expression, threshold = 0))
  } else {
    data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
      summarize(avg.exp = mean(log2(expm1(x = expression) + 1)), pct.exp = PercentAbove(x = expression, threshold = 0))
  }
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    mutate(avg.exp = MinMax(data = avg.exp, max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                    levels = sub(pattern = "-", replacement = ".", 
                                                 x = genes.plot), ordered = TRUE)
  
  #data.to.plot$pct.exp[data.to.plot$pct.exp <= dot.min] <- NA
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                 y = id)) + geom_point(mapping = aes(size = pct.exp, 
                                                                                     fill = avg.exp), colour="black",pch=21) + scale_radius(range = c(0, 
                                                                                                                                                      
                                                                                                                                                      
                                                                                                                                                      dot.scale), limits = c(0.25,1)) + scale_fill_gradientn(colors = cols.use, guide = guide_colorbar(ticks = FALSE)) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + coord_flip()
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, 
                                              hjust = 1))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}