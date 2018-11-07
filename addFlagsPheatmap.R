### Function to add flag

#### TO DO!!!!
#### STILL NEED TO CHECK OUT THE GROUPING ARGUMENTS. NGRP_V = 1 AND NGRP_V = 2 appear to be the same....

addFlag <- function(pheatmap,
                    labels_v,
                    repelDegree_v,
                    nGrp_v = 1) {
  #' Add special labels to pheatmap
  #' @description Add labels with flags to select rows of a heatmap
  #' @param pheatmap pheatmap object
  #' @param labels_v vector of labels to include
  #' @param repelDegree_v value between 0 and 1 that controls space to allocate for label repel.
  #' repel.degree = 0: spread out labels over existing range of kept labels
  #' repel.degree = 1: spread out labels of full y-axis
  #' @param nGrp_v numeric value to set the number of groups of flags. Flags within a group will be equidistant from each other
  #' @details Function provided by Z.Lin here https://stackoverflow.com/questions/52599180/partial-row-labels-heatmap-r
  #' @value pheatmap object
  #' @export
  
  ## Grab heatmap gtable
  heat_gtable <- pheatmap$gtable
  
  ## Extract row names grob
  label_grob <- heat_gtable$grobs[[which(heat_gtable$layout$name == "row_names")]]
  
  ## Replace all labels not found in labels_v with ''
  label_grob$label <- ifelse(label_grob$label %in% labels_v, label_grob$label, "")
  
  ## If there are more than 40 rows, space labels out by adding flags
  if (length(label_grob$label) > 40) {
    
    ## Calculate new y positions
    newY_v <- repelledY(d = label_grob$y,
                        d.select = label_grob$label != "",
                        k = repelDegree_v,
                        n = nGrp_v)
    
    ## Create new grob for the flag line segments
    newFlag_grob <- segmentsGrob(x0 = label_grob$x - unit(0.05, "npc"),
                                 x1 = label_grob$x + unit(0.25, "npc"),
                                 y0 = label_grob$y[label_grob$label != ""],
                                 y1 = newY_v)
    
    ## Shift selected labels to the right to make room for line segments
    label_grob$x <- label_grob$x + unit(0.3, "npc")
    
    ## Change y positions of selected labels
    label_grob$y[label_grob$label != ""] <- newY_v
    
    ## Add flag line segment grob to heatmap
    heat_gtable <- gtable::gtable_add_grob(x = heat_gtable,
                                           grobs = newFlag_grob,
                                           t = 4, l = 4)
  } # fi
  
  ## Replace label grob with
    ## nrow(heatmap) > 40: selected labels only, spaced evenly and with flags added
    ## nrow(heatmap) <= 40: selected labels only, with same y positioning
  heat_gtable$grobs[[which(heat_gtable$layout$name == "row_names")]] <- label_grob
  
  ## plot result
  grid.newpage()
  grid.draw(heat_gtable)
  
  ## Return copy invisibly
  invisible(heat_gtable)
}


### New repelledy Y
repelledY <- function(d, d.select, k = repelDegree_v, n = nGrp_v){
  #' Calculate y positioning
  #' @description Calculate evenly-spaced heatmap label y positions
  #' @param d vector of distances for labels (usually $y element of label_grob)
  #' @param d.select vector of T/F for which labels are significant (i.e. to be plotted)
  #' @param k repelDegree_v
  #' @param n nGrp_v
  #' @export
  
  ## Get full y-axis range of points
  fullYRange_v <- sapply(seq_along(d), function(i) as.numeric(convertX(d[i], "npc")))
  
  ## Subset for just the labeled ones
  selectYRange_v <- sapply(seq_along(d[d.select]), function(i) as.numeric(convertX(d[d.select][i], "npc")))
  
  ## Split select range
  chunks <- function(x, l) split(x, cut(seq_along(x), l, labels = F))
  if (n >= 2) {
    selectYRange_lsv <- chunks(selectYRange_v, n)
  } else {
    selectYRange_lsv <- list(selectYRange_v)
  }
  
  ## Get final values for each range - to do - make this work with repelDegree_v
  allFinal_lsv <- list()
  for (i in 1:length(selectYRange_lsv)) {
    currSelect_v <- selectYRange_lsv[[i]]
    currSmin_v <- min(currSelect_v); currSmax_v <- max(currSelect_v)
    # currFullMax_v <- max(fullYRange_v[fullYRange_v >= currSmax_v])
    # currFullMin_v <- min(fullYRange_v[fullYRange_v < currSmin_v])
    # currFinal_v <- unit(seq(from = currSmax_v + k * (currFullMax_v = currSmax_v),
    #                         to = currSmin_v - k * (currSmin_v - currFullMin_v),
    #                         length.out = length(currSelect_v)), "npc")
    currFinal_v <- unit(seq(from = currSmax_v, to = currSmin_v, length.out = length(currSelect_v)), "npc")
    allFinal_lsv[[i]] <- currFinal_v
  } # for i
  
  ## Change to npc
  allFinal_v <- unit(unlist(allFinal_lsv), "npc")
  
  ## Return
  return(allFinal_v)
} # repelledY














### Original repelledY Function with strip.npc
origRepelledY <- function(d, d.select, k = repelDegree_v){
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