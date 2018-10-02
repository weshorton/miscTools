### Function to add flag
addFlag <- function(pheatmap,
                    labels_v,
                    repelDegree_v) {
  #' Add special labels to pheatmap
  #' @description Add labels with flags to select rows of a heatmap
  #' @param pheatmap pheatmap object
  #' @param labels_v vector of labels to include
  #' @param repelDegree_v value between 0 and 1 that controls space to allocate for label repel.
  #' repel.degree = 0: spread out labels over existing range of kept labels
  #' repel.degree = 1: spread out labels of full y-axis
  #' @details Function provided by Z.Lin here https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
  #' @value pheatmap object
  #' @export
  
  ## Grab heatmap gtable
  heat_gtable <- pheatmap$gtable
  
  ## Extract row names grob
  label_grob <- heat_gtable$grobs[[which(heat_gtable$layout$name == "row_names")]]
  
  ## Replace all labels not found in labels_v with ''
  label_grob$label <- ifelse(label_grob$label %in% labels_v, label_grob$label, "")
  
  ## Calculate y positioning
  repelledY <- function(d, d.select, k = repelDegree_v){
    #' Calculate y positioning
    #' @description Calculate evenly-spaced heatmap label y positions
    #' @param d vector of distances for labels (usually $y element of label_grob)
    #' @param d.select vector of T/F for which labels are significant (i.e. to be plotted)
    #' @param k repelDegree_v
    #' @export
    
    ## Get current label positions (recursive)
    ## Unit is 'npc'
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      } # fi
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    } # strip.npc
    
    ## Get full y-axis range of points
    fullYRange_v <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    
    ## Subset for chosen k/repelDegree_v
    selectYRange_v <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    ## Get final values spaced by k
    finalVals_v <- unit(seq(from = max(selectYRange_v) + k*(max(fullYRange_v) - max(selectYRange_v)),
                            to = min(selectYRange_v) - k*(min(selectYRange_v) - min(fullYRange_v)),
                            length.out = sum(d.select)), "npc")
    ## Return
    return(finalVals_v)
  } # repelledY
  
  ## Calculate new y positions
  newY_v <- repelledY(d = label_grob$y,
                      d.select = label_grob$label != "")
  
  ## Create new grob for the flag line segments
  newFlag_grob <- segmentsGrob(x0 = label_grob$x,
                               x1 = label_grob$x + unit(0.15, "npc"),
                               y0 = label_grob$y[label_grob$label != ""],
                               y1 = newY_v)
  
  ## Shift selected labels to the right to make room for line segments
  label_grob$x <- label_grob$x + unit(0.2, "npc")
  
  ## Change y positions of selected labels
  label_grob$y[label_grob$label != ""] <- newY_v
  
  ## Add flag line segment grob to heatmap
  heat_gtable <- gtable::gtable_add_grob(x = heat_gtable,
                                         grobs = newFlag_grob,
                                         t = 4, l = 4)
  
  ## Replace label positions
  heat_gtable$grobs[[which(heat_gtable$layout$name == "row_names")]] <- label_grob
  
  ## plot result
  grid.newpage()
  grid.draw(heat_gtable)
  
  ## Return copy invisibly
  invisible(heat_gtable)
}