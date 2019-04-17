set.seed(123)
library(tidyr)

time = as.Date('2009-01-01') + 0:9

wiki_1 <- data.frame(
  W = sample(1:1000,10,replace = T),
  X = sample(1:100,10,replace = T),
  Y = sample(1:10,10,replace = T),
  Z = sample(1:10,10, replace = T)
)

wiki_2 <- data.frame(
  A = sample(500:1000,10,replace = T),
  B = sample(90:100,10,replace = T),
  C = sample(1:10,10,replace = T),
  D = sample(1:10,10,replace = T),
  E = sample(1:20,10,replace = T)
)


df1 <- wiki_1
df2 <- wiki_2
n_v <- 2
ref_v <- "Y"
cutoff_v <- 11

selectColsByMedian <- function(df1, df2, ref_v, n_v, cutoff_v) {
  #' Select Columns By Median
  #' @description Select any number of columns from a test data.frame whose median value is
  #' close to the median value of a specified column from a reference data.frame. "Close to"
  #' is determined as the absolute value of the difference in medians being less thant he specified cutoff.
  #' Outputs a new data.frame containing the reference data.frame's test column and all matching columns
  #' from the test data.frame
  #' @param df1 reference data.frame
  #' @param df2 test data.frame
  #' @param ref_v column from reference data.frame to test against
  #' @param n_v number of columns from df2 to select
  #' @param cutoff_v value to use to determine if test columns' medians are close enough
  #' @return data.frame with 1 column from df1 and matching columns from df2

  ## Get median of ref
  med_v <- median(df1[,ref_v], na.rm = T)
  
  ## Get other medians
  otherMed_v <- apply(wiki_2, 2, function(x) median(x, na.rm = T))
  
  ## Get differences
  medDiff_v <- sapply(otherMed_v, function(x) abs(med_v - x))
  
  ## Get whoever is within range (and order them)
  inRange_v <- sort(medDiff_v[medDiff_v < cutoff_v])
  inRangeCols_v <- names(inRange_v)
  
  ## Select random sample, if needed
  if (length(inRangeCols_v) > n_v){
    whichRandom_v <- sample(1:length(inRangeCols_v), size = n_v, replace = F)
  } else {
    whichRandom_v <- 1:length(inRangeCols_v)
  }
  finalCols_v <- inRangeCols_v[whichRandom_v]
  
  ## Final output
  out_df <- cbind(df1[,ref_v], df2[,finalCols_v])
  colnames(out_df) <- c(ref_v, finalCols_v)
  
  ## Return
  return(out_df)
} # selectColsByMedian

### 3 matching columns, select 2
match3pick2_df <- selectColsByMedian(df1 = wiki_1, df2 = wiki_2, ref_v = "Y", n_v = 2, cutoff_v = 5)
match3pick2_df2 <- selectColsByMedian(df1 = wiki_1, df2 = wiki_2, ref_v = "Y", n_v = 2, cutoff_v = 5)

### 2 matching columns, select 2
match2pick2_df <- selectColsByMedian(df1 = wiki_1, df2 = wiki_2, ref_v = "Y", n_v = 2, cutoff_v = 4)
