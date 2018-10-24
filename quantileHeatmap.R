###
### Heatmap Function
###

quantileHeat <- function(data_df, col_df = NA, row_df = NA, annCol_lsv = NA, heatCol_v = NA, 
                         rev_v = T, scale_v = "row", colClust_v = T, rowClust_v = T, sortClust_v = T, qc_v = T, outDir_v = NA, ...) {
  #' Custom Heatmap Function
  #' @description This heatmap function uses quartile scaling to distribute colors more evenly among measurements.
  #' @param data_df data.frame with counts to plot. Rows = gene/analyte; columns = patient/sample.
  #' @param col_df data.frame with metadata for column annotations. Must have the same column names as data_df
  #' @param row_df data.frame with metadata for row annotations. Must have the same row names as data_df
  #' @param annCol_lsv list with color distinctions for treatments. Each element name must be a column in col_df or row_df.
  #' @param heatCol_v Can be one of 3 things: (1) vector of color identifiers, (2) name of RColorBrewer palette, (3) name of other color fxn (e.g. inferno)
  #' @param rev_v logical. TRUE - reverse order of colors. FALSE - keep original order
  #' @param scale_v one of "row", "col", or NA. Determines whether or not to scale data, and in which direction.
  #' @param rowClust_v logical. TRUE - cluster rows (see sortClust_v for their output). FALSE - don't cluster.
  #' @param colClust_v logical. TRUE - cluster cols (see sortClust_v for their output). FALSE - don't cluster.
  #' @param sortClust_v logical. TRUE - output sorted dendrograms. FALSE - output default dendrograms
  #' @param qc_v logical. TRUE - output QC plots of color mapping. FALSE - only output the final heatmap.
  #' @param outDir_v path to output directory for heatmap and qc plots (if specified). If blank, will print to stdout.
  #' @param ... other arguments passed to pheatmap, such as title, cellwidth, etc.
  #' @value 
  #' @export
  
  ########################
  ### DEFINE FUNCTIONS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ########################
  
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm = T)
    breaks[!duplicated(breaks)]
  }
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  ####################
  ### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ####################
  
  require(RColorBrewer)
  require(pheatmap)
  require(viridis)
  
  ############################
  ### HANDLE ... ARGUMENTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############################
  
  extraParams_ls = list(...)
  
  ### For testing
  #extraParams_ls <- list(cellwidth = 10)
  
  ####################
  ### HANDLE COLOR ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ####################
  
  ### Possible RColorBrewer names
  possibleColors_v <- rownames(brewer.pal.info)
  
  ### Evalulate colors
  if (is.na(heatCol_v[1])) {
    heatColors_v <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
    altColors_v <- colorRampPalette(brewer.pal(11, "RdBu"))(99)
  } else if (heatCol_v[1] %in% possibleColors_v) {
    ## Get max number of colors
    n_v <- brewer.pal.info[rownames(brewer.pal.info) == heatCol_v[1], "maxcolors"]
    ## Make 
    heatColors_v <- colorRampPalette(brewer.pal(n_v, heatCol_v))(100);
    altColors_v <- colorRampPalette(brewer.pal(n_v, heatCol_v))(99) 
  } else if (length(heatCol_v) == 1) {
    heatColors_v <- do.call(heatCol_v, list(100))
    altColors_v <- do.call(heatCol_v, list(99))
  } else {
    heatColors_v <- heatCol_v
    altColors_v <- heatColors_v[-1]
  } # fi
  
  ### Reverse
  if (rev_v) {
    heatColors_v <- rev(heatColors_v)
    altColors_v <- rev(altColors_v)
  }
  
  ### Add names
  names(heatColors_v) <- paste0("A", 1:length(heatColors_v))
  names(altColors_v) <- paste0("A", 1:length(altColors_v))
  
  ##################
  ### SCALE DATA ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##################
  
#  scale_v <- ifelse(scale_v == "row", 1, ifelse(scale_v == "col", 2, NA))
  
  if (!is.na(scale_v)) {
    if (scale_v == "row") {
      data_df <- t(scale(t(data_df)))
    } else if (scale_v == "col") {
      data_df <- scale(data_df)
    } else {
      stop("Invalid argument for scale_v")
    }
    # mean_v <- apply(data_df, scale_v, mean, na.rm = T)
    # sd_v <- apply(data_df, scale_v, sd, na.rm = T)
    # data_df <- (data_df - mean_v) / sd_v
  } # fi
  
  ##########################
  ### INITIAL DENDROGRAM ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ##########################
  
  ### Cluster columns
  unsort_col_hclust <- hclust(dist(t(data_df)))
  sort_col_hclust <- sort_hclust(unsort_col_hclust)
  
  ### Cluster rows
  unsort_row_hclust <- hclust(dist(data_df))
  sort_row_hclust <- sort_hclust(unsort_row_hclust)
  
  ### List
  clust_ls <- list("Unsorted Columns" = unsort_col_hclust, "Sorted Columns" = sort_col_hclust,
                   "Unsorted Rows" = unsort_row_hclust, "Sorted Rows" = sort_row_hclust)
  
  ### Plot
  if (qc_v) {
    if (!is.na(outDir_v)) pdf(file = file.path(outDir_v, "dendrograms.pdf"))
    sapply(names(clust_ls), function(x) plot(clust_ls[[x]], main = x, xlab = "", sub = ""))
    if (!is.na(outDir_v)) graphics.off()
  } # fi
  
  #######################
  ### INITIAL HEATMAP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #######################
  
  if (qc_v) {
    ## Get file name
    fileName_v <- ifelse(is.na(outDir_v), NA, file.path(outDir_v, "initialHeat.pdf"))
    
    ## Get standard arguments
    stdArgs_v <- list(mat               = data_df,
                      cluster_rows      = unsort_row_hclust,
                      annotation_row    = row_df,
                      cluster_cols      = unsort_col_hclust,
                      annotation_col    = col_df,
                      annotation_colors = annCol_lsv,
                      color             = heatColors_v,
                      filename          = fileName_v)
    
    ## If want sorted cluster, add now
    if (sortClust_v) {
      stdArgs_v$cluster_rows <- sort_row_hclust
      stdArgs_v$cluster_cols <- sort_col_hclust
    } # fi
    
    ## If you don't want to cluster, add now
    if (!rowClust_v) stdArgs_v$cluster_rows <- F
    if (!colClust_v) stdArgs_v$cluster_cols <- F
    
    ## Plot with extra arguments
    do.call(pheatmap, c(stdArgs_v, extraParams_ls))
  } # fi
  
  ###########################
  ### TEST UNIFORM BREAKS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###########################
  
  ### Convert to matrix
  data_mat <- as.matrix(data_df)
  
  ### Plot density of values
  orig_density_df <- data.frame(values = as.numeric(data_mat))
  orig_density_gg <- ggplot(orig_density_df, aes(values)) + geom_density(bw = "SJ") + theme_bw(base_size = 16)
  
  ### Construct data to visualize color distribution
  orig_breaks_v <- seq(min(data_mat, na.rm = T), max(data_mat, na.rm = T), length.out = length(heatColors_v))
  orig_color_df <- data.frame(xmin_v = orig_breaks_v[1:(length(orig_breaks_v)-1)],
                              xmax_v = orig_breaks_v[2:length(orig_breaks_v)],
                              ymin_v = 0, ymax_v = max(density(data_mat, bw = "SJ", na.rm = T)$y),
                              fill_v = names(altColors_v),
                              stringsAsFactors = F)
  
  ### Plot density values with color distribution as well
  orig_color_density_gg <- ggplot() +
    geom_rect(data = orig_color_df, mapping = aes(xmin = xmin_v, xmax = xmax_v, ymin = ymin_v, ymax = ymax_v, fill = fill_v)) +
    geom_density(data = orig_density_df, mapping = aes(values), bw = "SJ", color = "cyan") +
    scale_fill_manual(values = altColors_v, breaks = names(altColors_v)) +
    labs(title = "Uniform Breaks") + my_theme + theme(legend.position = "none")
  
  ### Plot number of data points per color
  orig_pointsPer_df <- as.data.frame(table(cut(data_mat, orig_breaks_v)))
  orig_pointsPer_df$fill <- names(altColors_v)
  
  orig_pointsPer_gg <- ggplot() +
    geom_bar(data = orig_pointsPer_df, mapping = aes(x = Var1, weight = Freq, fill = fill), color = "black", size = 0.1) +
    coord_flip() +
    scale_fill_manual(values = altColors_v) +
    my_theme + theme(legend.position = "none", axis.text.y = element_text(size = 8, angle = 45)) +
    labs(y = "Data Points", x = "Breaks", title = "Uniform Breaks - Number of Data Points per Color")
  
  if (qc_v) {
    if (!is.na(outDir_v)) pdf(file = file.path(outDir_v, "uniformBreaks.pdf"))
    print(orig_density_gg)
    print(orig_color_density_gg)
    print(orig_pointsPer_gg)
    if (!is.na(outDir_v)) graphics.off()
  } # fi
  
  ############################
  ### TEST QUANTILE BREAKS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ############################
  
  ### Get new breaks
  breaks_v <- length(orig_breaks_v) + 1
  quantile_breaks_v <- quantile_breaks(data_mat, n = breaks_v)
  while(length(quantile_breaks_v) != breaks_v) {
    breaks_v <- breaks_v + 1
    quantile_breaks_v <- quantile_breaks(data_mat, n = breaks_v)
  } # while
  
  ### Construct color distribution data.frame
  quantile_color_df <- data.frame(xmin_v = quantile_breaks_v[1:(length(quantile_breaks_v)-1)],
                                  xmax_v = quantile_breaks_v[2:length(quantile_breaks_v)],
                                  ymin_v = 0, ymax_v = max(density(data_mat, bw = "SJ", na.rm = T)$y),
                                  fill_v = names(heatColors_v),
                                  stringsAsFactors = F)
  
  ## Plot density values with color distribution
  quantile_color_density_gg <- ggplot() +
    geom_rect(data = quantile_color_df, mapping = aes(xmin = xmin_v, xmax = xmax_v, ymin = ymin_v, ymax = ymax_v, fill = fill_v)) +
    geom_density(data = orig_density_df, mapping = aes(values), bw = "SJ", color = "cyan") +
    scale_fill_manual(values = heatColors_v, breaks = names(heatColors_v)) +
    labs(title = "Quantile Breaks") + my_theme + theme(legend.position = "none")
  
  ## Plot number of data points per color
  quantile_pointsPer_df <- as.data.frame(table(cut(data_mat, quantile_breaks_v)))
  quantile_pointsPer_df$fill <- names(heatColors_v)
  
  quantile_pointsPer_gg <- ggplot() +
    geom_bar(data = quantile_pointsPer_df, mapping = aes(x = Var1, weight = Freq, fill = fill), color = "black", size = 0.1) +
    coord_flip() +
    scale_fill_manual(values = heatColors_v) +
    my_theme + theme(legend.position = "none", axis.text.y = element_text(size = 8, angle = 45)) +
    labs(y = "Data Points", x = "Breaks", title = "Quantile Breaks - Number of Data Points per Color")
  
  if (qc_v) {
    if (!is.na(outDir_v)) pdf(file = file.path(outDir_v, "quantileBreaks.pdf"))
    print(quantile_color_density_gg)
    print(quantile_pointsPer_gg)
    if (!is.na(outDir_v)) graphics.off()
  } # fi
  
  ########################
  ### QUANTILE HEATMAP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ########################
  
  ## Get file name
  fileName_v <- ifelse(is.na(outDir_v), NA, file.path(outDir_v, "quantileHeat.pdf"))
    
  ## Get standard arguments
  stdArgs_v <- list(mat               = data_df,
                    cluster_rows      = unsort_row_hclust,
                    annotation_row    = row_df,
                    cluster_cols      = unsort_col_hclust,
                    annotation_col    = col_df,
                    annotation_colors = annCol_lsv,
                    color             = heatColors_v,
                    breaks            = quantile_breaks_v,
                    filename          = fileName_v)
  
  if (sortClust_v) {
    stdArgs_v[["cluster_rows"]] <- sort_row_hclust
    stdArgs_v[["cluster_cols"]] <- sort_col_hclust
  }
  
  ## If you don't want to cluster, add now
  if (!rowClust_v) stdArgs_v$cluster_rows <- F
  if (!colClust_v) stdArgs_v$cluster_cols <- F
  
  ## Plot with extra arguments
  do.call(pheatmap, c(stdArgs_v, extraParams_ls))

} # myHeat


# ###
# ### Test
# ###
# 
# ### Dependency
# source("~/my_tool_repos/WesPersonal/utilityFxns.R")
# 
# ### Data
# data_df <- convertDFT(fread("~/projs/Schedin/HMIBU/input/cluster_CD45-HM-NSAID_082018.txt"))
# data_df <- data_df[!(rownames(data_df) %in% c("Intratumoral CD45+", "Tumor_Border CD45+")),]
# 
# ### Meta
# col_df <- convertDFT(fread("~/projs/Schedin/HMIBU/input/meta.txt"))
# row_df <- data.frame("Site" = gsub(" .*", "", rownames(data_df)))
# rownames(row_df) <- rownames(data_df)
# 
# ### Colors
# colors_v <- brewer.pal(11, "Spectral")
# groupColors_v <- colors_v[9:11]; names(groupColors_v) <- c("INV", "INVIBU", "NP")
# siteColors_v <- colors_v[c(1,3,5)]; names(siteColors_v) <- c("Intratumoral", "Tumor_Border", "Unknown")
# annCol_lsv <- list("Group" = groupColors_v,
#                    "Site" = siteColors_v)
# 
# ### Other args
# heatCol_v <- "RdBu"
# rev_v <- T
# scale_v <- "row"
# sortClust_v <- T
# qc_v <- T
# outDir_v <- "~/Desktop/heatTest/"
# 
# quantileHeat(data_df = data_df,
#        col_df = col_df,
#        row_df = row_df,
#        annCol_lsv = annCol_lsv,
#        heatCol_v = heatCol_v,
#        rev_v = rev_v,
#        scale_v = scale_v,
#        sortClust_v = sortClust_v,
#        qc_v = qc_v,
#        outDir_v = outDir_v,
#        cellwidth = 10,
#        main = "Test Heatmap")