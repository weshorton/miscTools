### From https://github.com/kevinblighe/EnhancedVolcano/blob/master/R/EnhancedVolcano.R

### Uncomment this to assign all the standard arguments
# xlim = c(min(toptable[,x], na.rm=TRUE),
#          max(toptable[,x], na.rm=TRUE))
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


EnhancedVolcano <- function(
  toptable,
  lab,
  x,
  y,
  selectLab = NULL,
  plotLabels = NA,
  both = T, # plot both significant labels and selected labels
  colorSelect = T, # change the color of select label **points** to red4
  xlim = c(min(toptable[,x], na.rm=TRUE),
           max(toptable[,x], na.rm=TRUE)),
  ylim = c(0, max(-log10(toptable[,y]), na.rm=TRUE) + 5),
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~italic(P)),
  axisLabSize = 16,
  pCutoff = 0.05,
  pLabellingCutoff = pCutoff,
  FCcutoff = 2.0,
  title = "",
  titleLabSize = 16,
  transcriptPointSize = 0.8,
  transcriptLabSize = 2.0,
  col = c("grey30", "forestgreen", "royalblue", "red2"),
  colAlpha = 1/2,
  legend = c("NS","Log2 FC","P","P & Log2 FC"),
  legendPosition = "top",
  legendLabSize = 10,
  legendIconSize = 3.0,
  DrawConnectors = FALSE,
  widthConnectors = 0.5,
  colConnectors = "black",
  cutoffLineType = "longdash",
  cutoffLineCol = "black",
  cutoffLineWidth = 0.4)
{
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
  
  ### Create significance groupings
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[,x]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[,y]<pCutoff)] <- "P"
  toptable$Sig[(toptable[,y]<pCutoff) &
                 (abs(toptable[,x])>FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig,
                         levels=c("NS","FC","P","FC_P"))
  
  ### Add label column and plotting values
  if (length(lab) == nrow(toptable)){
    toptable$lab <- lab
  } else {
    ## Get the ones that are specified in labs
    labs <- sapply(rownames(toptable), function(x) {
      y <- ifelse(x %in% lab, x, "")
    })
    toptable$lab <- labs
    ## Also need to get the ones specified in selectLab
    for (i in 1:length(selectLab)){
      toptable[rownames(toptable) == selectLab[i], "lab"] <- selectLab[i]
    }
  }
  toptable$xvals <- toptable[,x]
  toptable$yvals <- toptable[,y]
  
  # ### Add separate color for specified labels
  # ### Also turn significance groupings into factor
  # if (colorSelect) {
  #   toptable$Sig[toptable$lab %in% selectLab] <- "S"
  #   # toptable$Sig <- factor(toptable$Sig,
  #   #                        levels = c("NS", "FC", "P", "FC_P", "S"))
  #   toptable$Sig <- factor(toptable$Sig,
  #                          levels = c("S", "NS", "FC", "P", "FC_P"))
  #   #col <- c(col, "red4")
  #   #legend <- c(legend, "Select")
  #   col <- c("red4", col)
  #   legend <- c("Select", legend)
  # } else {
  #   toptable$Sig <- factor(toptable$Sig,
  #                          levels=c("NS","FC","P","FC_P"))
  # }
  
  ### Add names to col and legend
  names(col) <- levels(toptable$Sig)
  names(legend) <- levels(toptable$Sig)
  
  if (min(toptable[,y], na.rm=TRUE) == 0) {
    warning("One or more P values is 0. Converting to minimum possible value...", call. = FALSE)
    toptable[which(toptable[,y] == 0), y] <- .Machine$double.xmin
  }
  
  ### Replace labels with selected labels, if specified.
  if (!is.null(selectLab) & !both) {
    names.new <- rep("", length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  
  ## Subset data, only used if !is.null(selectLab)
  subTopTable <- subset(toptable, toptable[,y] < pLabellingCutoff & abs(toptable[,x]) > FCcutoff)
  subTopTable <- rbind(subTopTable, toptable[toptable$lab %in% selectLab,])
  subTopTable <- subTopTable[subTopTable$lab != "",]
  ## Add colors
  subTopTable$Color <- "black"
  subTopTable$Color[subTopTable$lab %in% selectLab] <- "#8B0000"

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
    plot <- plot + 
      annotate("text", label = plotLabels[1], x = xlim[1], y = ylim[1], fontface = "bold") +
      annotate("text", label = plotLabels[2], x = xlim[2], y = ylim[1], fontface = "bold")
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
    ## Subset data
    subTopTable <- subset(toptable, toptable[,y] < pLabellingCutoff & abs(toptable[,x])>FCcutoff)
    subTopTable <- rbind(subTopTable, toptable[toptable$lab %in% selectLab,])
    subTopTable <- subTopTable[subTopTable$lab != "",]
    ## Add colors
    subTopTable$Color <- "black"
    subTopTable$Color[subTopTable$lab %in% selectLab] <- "#8B0000"
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