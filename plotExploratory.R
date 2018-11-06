plotExploratory <- function(dir_v,            # Path to output directory
                            name_v,           # Name for file and titles
                            data,             # object for plotting (SeqExpressionSet)
                            colors_lsv,       # List of colors, one element must be treat_v
                            meta_df,          # metadata, one column must be treat_v
                            treat_v,          # Treatment distinction for colors, plotting, etc.
                            shape_v = NULL,   # Secondary treatment/division for PCA points
                            type_v,           # Character vector - Raw or Norm
                            labels_v = T,     # logical. TRUE - use colnames as labels of PCA points. FALSE - just use points
                            plot_v            # logical. TRUE - write to file. FALSE - write to console
                            ) {
  
  ### Note that to actually access the counts of a normalized SeqExpressionSet, you have to use normCounts(data)
  ### counts(rawData) and counts(normData) are the same, but normCounts(normData) is different.
  ### The above does not matter for this function, however, because EDASeq::plotRLE and EDASeq::plotPCA
  ### recognize a normalized SeqExpressionSet and use the normalized counts.
  ### Get color vector
  colors_v <- colors_lsv[[treat_v]][meta_df[[treat_v]]]
  
  ### Get shape vector
  if (!is.null(shape_v)) {
    ### Get different pch values
    shapes_v <- c(20,17,15,18,1,2,0,5)
    ### Get table of occurrences
    nShapes_v <- as.data.table(table(meta_df[[shape_v]]))
    ### Subset shapes and name them
    shapes_v <- shapes_v[1:nrow(nShapes_v)]; names(shapes_v) <- nShapes_v$V1
    ### Map
    pch_v <- shapes_v[match(meta_df[[shape_v]], names(shapes_v))]
  } else {  
  ### Plot boxplots
  if (plot_v) pdf(file = file.path(dir_v, paste(name_v, type_v, "RLE.pdf", sep = "_")))
  plotRLE(data, outline = F, col = colors_v,
          main = paste0("RLE: Relative Log Expression\n", type_v, " Counts - ", name_v),
          xlab = "Sample", ylab = "log(Read Count / Median Count Across Samples)", xaxt = 'n')
  legend("bottomright", legend = levels(meta_df[[treat_v]]), col = colors_lsv[[treat_v]], lwd = 2, bty = "n")
  Sys.sleep(1)
  if (plot_v) dev.off()
  Sys.sleep(1)
  
  ### Plot PCA
  if (plot_v) pdf(file = file.path(dir_v, paste(name_v, type_v, "PCA.pdf", sep = "_")))
  par(mar = c(5.1, 5.1, 4.1, 8.1))
  par(xpd = T)
  EDASeq::plotPCA(data, col = colors_v, labels = labels_v, pch = pch_v, main = paste0("PCA of ", type_v, " Counts - ", name_v))
  legend("topright", inset = c(-0.26, 0), legend = levels(meta_df[[treat_v]]), col = colors_lsv[[treat_v]], 
         pch = 20, lty = NA, lwd = 2, bty = 'n')
  if (!is.null(shape_v)) legend("right", legend = names(shapes_v), inset = c(-0.09, 0), pch = shapes_v, bty = 'n')
  Sys.sleep(1)
  if (plot_v) dev.off()
  Sys.sleep(1)
  par(xpd = F)
}
