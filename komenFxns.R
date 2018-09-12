###
### KOMEN FUNCTIONS
###

######## Table of Contents ########
### splitInput                  ###
### ratioInput                  ###
###################################

### Functions commonly used for Komen Analysis

##################
### Split Data ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################

### Split the combined input table into a list with one table per measurement
splitInput <- function(combined_dt, split_v, metaCols_v, measureCols_v){
  #' Split a large data.table into a list of multiple subset data.frames
  #' @description Split data.table by measurement time points into list of subset data.frames
  #' @param combined_dt Large data.table with non-unique patient rows (i.e. each patient has multiple rows, for each time point)
  #' @param split_v time point column name. Where the identities of the multiple time points reside. Used to subset data.
  #' @param metaCols_v metadata columns that should be ignored while subsetting
  #' @param measureCols_v measurement columns that need to be included
  #' @value Export a list of data.frames containing the non-NA measurements of each division specified in split_v.
  #' @export
  
  ## NOTE (2018-06-15) I'M NOT SURE WHAT MEASURECOLS_V IS DOING. THESE COLUMNS (E.G. OUTCOMETYPE) DON'T
  ## EXIST IN THE INPUT DATA. THEY'RE ONLY IN THE METADATA.
  
  ## Get divisions
  divisions_v <- unique(combined_dt[[split_v]])
  
  ## Split for each
  splitData_lsdf <- sapply(divisions_v, function(x) {
    ## Subset data
    currSub_dt <- combined_dt[get(split_v) == x,]
    ## Remove meta info and change to data.frame
    currExpr_df <- as.data.frame(currSub_dt[,!which(colnames(currSub_dt) %in% metaCols_v), with = F])
    ## Re-add row.names
    rownames(currExpr_df) <- currSub_dt[[1]]
    ## Add measure column info
    for (col_v in measureCols_v) currExpr_df[[col_v]] <- currSub_dt[[col_v]]
    ## Export
    return(currExpr_df)
  }, simplify = F)
} # splitInput

### Convert data into ratios
ratioInput <- function(splitData_lsdf, split_v, baseline_v = "Screen", measureCols_v = "foo") {
  #' Turn data into ratios
  #' @description Turn all non-baseline values into a ratio b/w that timepoint and baseline (formula log2(timepoint / baseline))
  #' @param splitData_lsdf list of data.frames. One data.frame for each measurement point. Output of splitInput()
  #' @param split_v Values from time point column. Should be a vector of unique time points. Doesn't have to include baseline, but can.
  #' @param measureCols_v measurement columns that need to be included
  #' @param baseline_v Element of split_v. Baseline measurement. Default is "Screen".
  #' @value list of data.frames that is one element shorter than splitData_lsdf. Values are now ratios instead of counts.
  #' @export
  
  ## NOTE (2018-06-15) I'M NOT SURE WHAT MEASURECOLS_V IS DOING. THESE COLUMNS (E.G. OUTCOMETYPE) DON'T
  ## EXIST IN THE INPUT DATA. THEY'RE ONLY IN THE METADATA.
  
  ## Get divisions and remove baseline
  divisions_v <- grep(baseline_v, split_v, value = T, invert = T)
  
  ## Split data
  baseline_df <- splitData_lsdf[[baseline_v]]
  values_lsdf <- splitData_lsdf[divisions_v]
  
  ## Turn value_lsdf to ratios
  ratios_lsdf <- sapply(divisions_v, function(x) {
    
    cat(sprintf("Currently on division: %s\n", x))
    
    ## Get data
    currData_df <- values_lsdf[[x]]
    
    ## Find patients in baseline, but not current value data.frame
    rmBaseline_v <- paste(setdiff(rownames(baseline_df), rownames(currData_df)), collapse = "|")
    
    ## Find patients in value data.frame, but not in baseline
    rmValue_v <- paste(setdiff(rownames(currData_df), rownames(baseline_df)), collapse = "|")
    
    ## Remove them
    if (rmValue_v != ""){
      currData_df <- currData_df[grep(rmValue_v, rownames(currData_df), invert = T, value = T),]
      cat(sprintf("Removed %s patient(s) from %s because not in %s\n\n", rmValue_v, x, baseline_v))
    }
    if (rmBaseline_v != ""){
      currBaseline_df <- baseline_df[grep(rmBaseline_v, rownames(baseline_df), invert = T, value = T),]
      cat(sprintf("Removed %s patient(s) from %s because not in %s\n\n", rmBaseline_v, baseline_v, x))
    }
    
    ## Get non-measurement columns (i.e. columns to take ratios of)
    currRatioCols_v <- grep(paste(measureCols_v, collapse = "|"), colnames(currData_df), invert = T, value = T)
    
    ## Turn into Ratios
    currRatio_df <- log2( as.matrix(currData_df[,currRatioCols_v]) / as.matrix(currBaseline_df[,currRatioCols_v]) )
    
    ## Change back to data.frame and add back measurement columns
    currRatio_df <- as.data.frame(currRatio_df)
    for (col_v in measureCols_v) currRatio_df[[col_v]] <- currData_df[[col_v]]
    
    ## Export
    return(currRatio_df)
  }, simplify = F)
} # ratioInput

### Get top 10 combinations that have a significant survival analysis p-value
getTop10 <- function(survResults_dt, comboCols_v, nameCols_v, pthresh_v){
  #' Subset for top 10 survival combinations
  #' @description Given a data.table of gene combinations and their survival analysis p-value results, output the
  #' top 10 most-significant combinations.
  #' @param survResults_dt data.table containing one row for each gene combination. First 1-3 columns are the gene names for the combo,
  #' second set of columns are the p-value results for the various +/- combinations. Final column is plot name
  #' @param comboCols_v column names for the different +/- combinations of the genes in each row
  #' @param nameCols_v column names for the columns containing the genes as well as the plot names
  #' @param pthresh_v p-value threshold to determine significance
  #' @value output a list of two data.tables. The first has all of the significant results for a given column. The second
  #'  has the same columns as the input, but subset for only the 10 most significant rows.
  #' @export
  
  ## Extract all significant results for each column individually
  individualSigResults_lsdt <- sapply(comboCols_v, function(x) {
    currSigData_dt <- survResults_dt[ get(x) < pthresh_v, mget(c(nameCols_v, x))]}, simplify = F)
  
  ## Subset each of those data.tables for the top10 only
  individualTop10Results_lsdt <- sapply(comboCols_v, function(y) {
    ## Get data
    currData_dt <- individualSigResults_lsdt[[y]]
    ## Subset
    if (nrow(currData_dt) > 10){
      ## Rank p-values
      currRank_v <- rank(currData_dt[[y]])
      ## Get indexes of top 10
      currRankIndex_v <- which(currRank_v <= 10)
      counter_v <- 1
      ## I don't think this does anything
      # while(length(currRankIndex_v) < 10){
      #   currRankIndex_v <- which(currRank_v < (10+counter_v))
      #   counter_v <- counter_v + 1
      # }
      ## Make data.dable
      currTop_dt <- currData_dt[currRankIndex_v,]
    } else {
      currTop_dt <- currData_dt
    } # fi
    ## Sort
    currTop_dt <- currTop_dt[order(get(y))]
    ## Export
    return(currTop_dt)
  }, simplify = F)
  
  ## Return
  return(list("allResults" = individualSigResults_lsdt, "top10Results" = individualTop10Results_lsdt))
  
}