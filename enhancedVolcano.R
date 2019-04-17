### From https://github.com/kevinblighe/EnhancedVolcano/blob/master/R/EnhancedVolcano.R

### Uncomment this to assign all the standard arguments
# xlim = c(min(toptable[,x], na.rm=TRUE),
#          max(toptable[,x], na.rm=TRUE))
# ylim = c(0, max(-log10(toptable[,y]), na.rm=TRUE) + 5)
# xlab = bquote(~Log[2]~ "fold change")
# axisLabSize = 16
# pLabellingCutoff = pCutoff
# titleLabSize = 16
# col = c("grey30", "forestgreen", "royalblue", "red2")
# colAlpha = 0.5
# legend = c("NS","Log2 FC","P","P & Log2 FC")
# legendLabSize = 10
# legendIconSize = 3.0
# widthConnectors = 0.5
# colConnectors = "black"
# cutoffLineType = "longdash"
# cutoffLineCol = "black"
# cutoffLineWidth = 0.4
# selectLab = NULL # watch this - often assigned in code
# both = T
# colorSelect = F


EnhancedVolcano <- function(toptable, lab, labCol = NULL, x, y, selectLab = NULL, 
                            plotLabels = NA, both = T, colorSelect = T,
                            xlim = c(min(toptable[,x], na.rm=TRUE), max(toptable[,x], na.rm=TRUE)),
                            ylim = c(0, max(-log10(toptable[,y]), na.rm=TRUE) + 5),
                            xlab = bquote(~Log[2]~ "fold change"), ylab = bquote(~-Log[10]~italic(P)),
                            axisLabSize = 16, pCutoff = 0.05, pLabellingCutoff = pCutoff, FCcutoff = 2.0, 
                            title = "", titleLabSize = 16, transcriptPointSize = 0.8, transcriptLabSize = 2.0, 
                            col = c("grey30", "forestgreen", "royalblue", "red2"), colAlpha = 1/2,
                            legend = c("NS","Log2 FC","P","P & Log2 FC"), legendPosition = "top", legendLabSize = 10, 
                            legendIconSize = 3.0, DrawConnectors = FALSE, widthConnectors = 0.5, colConnectors = "black", 
                            cutoffLineType = "longdash", cutoffLineCol = "black", cutoffLineWidth = 0.4) {
  #' Enhanced Volcano
  #' @description Plot pretty volcano plots using ggplot2. Allows for significant customization
  #' @param toptable some sort of table with results to be plot. Most common is a DESeqResults object.
  #' Rownames are genes, columns are results such as fold change, mean expression, etc.
  #' @param lab Character vector of genes to label. Can either be the rownames of toptable, or a subset of genes to label 
  #' (but still must be contained within rownames(toptable). CHECK THIS -> will only be labeld if above lfc and pvalue cut-off.
  #' @param labCol Character vector of column name to use to assign labels, instead of rownames. Default is NULL, which causes rownames to be used.
  #' if 'lab' is a subset of genes to label, they must be valid entries 
  #' @param x Column name from toptable that will be used for the x-axis values. Usually 'log2FoldChange'
  #' @param y Column name from toptable that will be used for the y-axis values. Usually 'pvalue' or 'padj'
  #' @param selectLab Vector of specific lables to call out with a different color. Will be labelled regardless of position on graph.
  #' @param plotLabels Vector of length two, specifying treatment labels for bottom left/right of plot. 
  #' First element is left-side lable, second is right-side
  #' @param both boolean determining if significant labels (i.e. above lfc and pvalue cut-off) AND selected labels should be plot (TRUE).
  #' or if just selected lables should be plot (FALSE). NEED TO IMPROVE THIS.
  #' @param colorSelect boolean determining if selected labels should be colored differently or not. Default TRUE, which colors select points to red4
  #' @param xlim Numeric vector of length two determining the x-axis limits
  #' @param ylim Numeric vector of length two determining the y-axis limits
  #' @param xlab Character vector for x-axis label
  #' @param ylab Character vector for y-axis label
  #' @param axisLabSize Numeric vector determining axis label size. Default 16.
  #' @param pCutoff Numeric vector determining the p-value (y-axis) value at which to draw dotted line
  #' @param pLabellingCutoff Numeric vector determining the p-value (y-axis) value at which to begin labelling significant points
  #' @param FCcutoff Numeric vector determining the lfc (x-axis) value at which to draw dotted line
  #' @param title Character vector specifying plot title
  #' @param titleLabSize Numeric vector specifying the title size. 
  #' @param transcriptPointSize Numeric vector determining cex of points. Default 0.8
  #' @param transcriptLabSize Numeric vector determining label size. Default 2
  #' @param col colors to use for points. Length of four. Order is: (1) below pvalue and lfc cut-off; (2) below pvalue and above lfc;
  #' (3) above pvalue and below lfc; (4) above pvalue and above lfc
  #' @param colAlpha level of transparency of points. Default 0.5
  #' @param legend Legend names for four divisions of points. See 'col' for description. Must be in this order.
  #' @param legendPosition Where to draw legend
  #' @param legendLabSize size of legend labels
  #' @param legendIconSize cex of legend points
  #' @param DrawConnectors boolean. Should connectors be drawn from points to labels?
  #' @param widthConnectors lwd of connector lines. Default 0.5
  #' @param colConnectors color of connector lines. Default 'black'
  #' @param cutoffLineType lty of pvalue and lfc cut-off lines. Default "longdash"
  #' @param cutoffLineCol color of cut-off lines. Default 'black'
  #' @param cutoffLineWidth lwd of cut-off lines. Default 0.4
  #' @details
  #' ## When to change 'lab'
  #' If there are many genes that are above the lfc and p-value cut-offs, the plot gets pretty confusing. Instead, you use more
  #' stringent cut-offs prior to plotting in order to determine a set of genes that are extremely significant that you would like to lable.
  #' 
  #' @export
  
  ###
  ### Check Dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  if(!requireNamespace("ggplot2")) {
    stop("Please install ggplot2 first.", call.=FALSE)
  }
  
  if(!requireNamespace("ggrepel")) {
    stop("Please install ggrepel first.", call.=FALSE)
  }
  
  if(!is.numeric(toptable[,x])) {
    stop(paste(x[i], " is not numeric!", sep=""))
  }
  
  if(!is.numeric(toptable[,y])) {
    stop(paste(x[i], " is not numeric!", sep=""))
  }
  
  requireNamespace("ggplot2")
  requireNamespace("ggrepel")
  
  ### Initialize variables
  i <- xvals <- yvals <- Sig <- NULL
  
  ### Create significance groupings using provided cut-offs
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[,x]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[,y]<pCutoff)] <- "P"
  toptable$Sig[(toptable[,y]<pCutoff) &
                 (abs(toptable[,x])>FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig,
                         levels=c("NS","FC","P","FC_P"))
  
  ###
  ### Add label column and plotting values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ## Simply re-assign to toptable, if all genes provided
  if (length(lab) == nrow(toptable)){
    toptable$lab <- lab
    
  ## If only a subset of genes are provided, have to put blanks in the others.
  ## Also need to include 'selectLab', if specified.
  ## Have to check if rownames or other column
  } else {
    
    ## If rownames
    if (is.null(labCol)) {
      
      ## Replace non-selected with ''
      labs <- sapply(rownames(toptable), function(x) {
        y <- ifelse(x %in% lab, x, "")
        })
      toptable$lab <- labs
      
      ## Also need to get the ones specified in selectLab
      for (i in 1:length(selectLab)) toptable[rownames(toptable) == selectLab[i], "lab"] <- selectLab[i]
      
    ## If selected column
    } else {
      
      ## New column  
      toptable$lab <- toptable[[labCol]]
      
      ## Just use indexes, if provided
      if (is.numeric(lab)) {
        
        toptable[-lab, "lab"] <- ""
        for (i in 1:length(selectLab)) toptable[selectLab[i], "lab"] <- toptable[selectLab[i],][[labCol]]
        
      ## Otherwise follow same as above, but column instead of rownames
      } else {
        
        labs <- sapply(toptable[[labCol]], function(x) {
          y <- ifelse(x %in% lab, x, "")
        })
        toptable$lab <- labs
        
        for (i in 1:length(selectLab)) toptable[toptable[[labCol]] == selectLab[i], "lab"] <- selectLab[i]
      } # fi
      
    } # fi
    
  } # fi
  
  ## New columns for x and y values so that column name is standardized
  toptable$xvals <- toptable[,x]
  toptable$yvals <- toptable[,y]
  
  ## Add names to col and legend
  names(col) <- levels(toptable$Sig)
  names(legend) <- levels(toptable$Sig)
  
  ## Handle 0 p-values
  if (min(toptable[,y], na.rm=TRUE) == 0) {
    warning("One or more P values is 0. Converting to minimum possible value...", call. = FALSE)
    toptable[which(toptable[,y] == 0), y] <- .Machine$double.xmin
  }
  
  ### If both is FALSE and selectLab is provided, just plot selectLabel. Make all other labels ''
  if (!is.null(selectLab) & !both) {
    names.new <- rep("", length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  
  ###
  ### Handle Subset Labels ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ## Subset data based on provided cut-offs
  subTopTable <- subset(toptable, toptable[,y] < pLabellingCutoff & abs(toptable[,x]) > FCcutoff)
  
  ## Find selectLab's that aren't already in subTopTable
  notSigSelectLab <- selectLab[!which(selectLab %in% rownames(subTopTable))]
  
  ## Add them to the table, if needed
  if (length(notSigSelectLab) > 0) {
    subTopTable <- rbind(subTopTable, toptable[toptable$lab %in% notSigSelectLab,])
  } # fi
  
  ## Remove empty labels. Labels would be empty if they passed the lfc and p-value cut-offs, but aren't in 'lab' argument
  subTopTable <- subTopTable[subTopTable$lab != "",]
  
  ## Add colors (only works if there is something there!)
  if (nrow(subTopTable) > 0) {
    subTopTable$Color <- "black"
    subTopTable$Color[subTopTable$lab %in% selectLab] <- "#8B0000"
  } else {
    subTopTable$Color <- character()
  }

  ## Begin plot
  plot <- ggplot2::ggplot(toptable,
                          ggplot2::aes(x=xvals, y=-log10(yvals))) +
    
    ggplot2::geom_point(ggplot2::aes(color=factor(Sig)),
                        alpha=colAlpha, size=transcriptPointSize) +
    
    ggplot2::scale_color_manual(values = col, labels = legend) + 
    
    ggplot2::theme_bw(base_size=24) +
    
    ggplot2::theme(
      legend.background=ggplot2::element_rect(),
      plot.title=ggplot2::element_text(angle=0,
                                       size=titleLabSize,
                                       face="bold",
                                       vjust=1,
                                       hjust =0.5),
      
      panel.grid.major=ggplot2::element_blank(),
      panel.grid.minor=ggplot2::element_blank(),
      
      axis.text.x=ggplot2::element_text(angle=0,
                                        size=axisLabSize,
                                        vjust=1),
      axis.text.y=ggplot2::element_text(angle=0,
                                        size=axisLabSize,
                                        vjust=1),
      axis.title=ggplot2::element_text(size=axisLabSize),
      
      legend.position=legendPosition,
      legend.key=ggplot2::element_blank(),
      legend.key.size=ggplot2::unit(0.5, "cm"),
      legend.text=ggplot2::element_text(
        size=legendLabSize),
      title=ggplot2::element_text(
        size=legendLabSize),
      legend.title=ggplot2::element_blank()) +
    
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes=list(size=legendIconSize))) +
    
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    
    ggplot2::xlim(xlim[1], xlim[2]) +
    ggplot2::ylim(ylim[1], ylim[2]) +
    
    ggplot2::ggtitle(title) +
    
    ggplot2::geom_vline(xintercept=c(-FCcutoff, FCcutoff),
                        linetype=cutoffLineType,
                        colour=cutoffLineCol,
                        size=cutoffLineWidth) +
    
    ggplot2::geom_hline(yintercept=-log10(pCutoff),
                        linetype=cutoffLineType,
                        colour=cutoffLineCol,
                        size=cutoffLineWidth)
  
  ## Add labels for each type
  if (!is.na(plotLabels[1])) {
    
    ## Determine limits
    lenPlotLabels_v <- nchar(plotLabels)
    limFactors_v <- sapply(lenPlotLabels_v, function(x) ifelse(x > 10, 0.85, 0.98))
    
    ## Add
    plot <- plot + 
      annotate("text", label = plotLabels[1], x = xlim[1]*limFactors_v[1], y = ylim[1], fontface = "bold") +
      annotate("text", label = plotLabels[2], x = xlim[2]*limFactors_v[2], y = ylim[1], fontface = "bold")
  }
  
  ## Add color for specified gene(s)
  if (colorSelect) {
    subSubTopTable <- subTopTable[subTopTable$Color == "#8B0000",]
    plot <- plot + ggplot2::geom_point(data = subSubTopTable,
                               ggplot2::aes(x = xvals, y = -log10(yvals)),
                               color = "red4")
  }
  
  if (DrawConnectors == TRUE && is.null(selectLab)) {
    plot <- plot + ggrepel::geom_text_repel(
      data=subset(toptable,
                  toptable[,y]<pLabellingCutoff &
                    abs(toptable[,x])>FCcutoff),
      ggplot2::aes(label=subset(toptable,
                                toptable[,y]<pLabellingCutoff &
                                  abs(toptable[,x])>FCcutoff)[,"lab"]),
      size = transcriptLabSize,
      segment.color = colConnectors,
      segment.size = widthConnectors,
      vjust = 1.0)
  } else if (DrawConnectors == TRUE && !is.null(selectLab)) {
    ## Add to plot
    plot <- plot + ggrepel::geom_text_repel(
      data = subTopTable,
      ggplot2::aes(label = lab),
      color = subTopTable$Color,
      size = transcriptLabSize,
      segment.colour = colConnectors,
      segment.size = widthConnectors,
      vjust = 1.0) 
  } else if (DrawConnectors == FALSE && !is.null(selectLab)) {
    ## Add to plot
    plot <- plot + ggplot2::geom_text(data = subTopTable,
                                      ggplot2::aes(label = lab),
                                      color = subTopTable$Color,
                                      size = transcriptLabSize,
                                      check_overlap = T, vjust = 1)
  } else if (DrawConnectors == FALSE && is.null(selectLab)) {
    plot <- plot + ggplot2::geom_text(data=subset(toptable,
                                                  toptable[,y]<pLabellingCutoff &
                                                    abs(toptable[,x])>FCcutoff),
                                      ggplot2::aes(label=subset(toptable,
                                                                toptable[,y]<pLabellingCutoff &
                                                                  abs(toptable[,x])>FCcutoff)[,"lab"]),
                                      size = transcriptLabSize,
                                      check_overlap = TRUE,
                                      vjust = 1.0)
  }
  
  return(plot)
}