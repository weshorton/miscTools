###
### Survival Functions
###

######## Table of Contents ########
### comboSurvFit                ###
### addSubCols                  ###
###################################

### Functions that help during survival analysis
comboSurvFit <- function(surv_obj, genes_v, data_df, cols_v = c("Days", "Response"), median_log = F, zScore_v = F, pthresh_v = 0.05, plot_v = T){
  #' Run Kaplan-Meier survival estimate on specified gene(s) expression data
  #' @description Given a base surival object, input data (expression fold-change from baseline), and a 1-3 gene set (columns from expression input), 
  #' run survival analysis on different combinations of those genes. Does so by adding/subtracting the gene sets in different configurations
  #' and then converting the expression data to 0 and 1. 0 if value below median, 1 if above.
  #' If p-value is below threshold, create survival plot. Returns p-values for all combinations.
  #' @param surv_obj object created via Surv(time,event) call
  #' @param genes_v vector of 2 or 3 genes to combine counts
  #' @param data_df count data where rows = patient/sample and columns = genes
  #' @param cols_v clinical information columns contained in data_df
  #' @param median_log TRUE - replace NA values of each gene with its median value. FALSE - keep NAs
  #' @param zScore_v TRUE - divide patients into 3 groups ("sigLow" (z < -1.96), "sigHigh" (z > 1.96), and "noChange" (-1.96 < z < 1.96)), after
  #' converting ratio data to z-scores. FALSE - use ratio data to divide groups into "low" (below median value) and "high" (above median value)
  #' @param pthresh_v p-value threshold to determine which survival plots to make
  #' @param plot_v TRUE - print plot; FALSE - no plotting
  #' @value Outputs list with 3 elements. (1) matrix containing plotting data (hi/low divisions and clinical info),
  #' (2) p-values of the different combinations, (3) calc-matrix to check hi/low designations. nrow(1) == length(2) == nrow(3) and
  #' also each index is the same.
  #' @export

  ## Subset data
  clin_df <- data_df[,cols_v]
  sub_df <- data_df[,genes_v, drop = F]
  
  ## Add median values, if specified
  if (median_log) {
    sub_df <- apply(sub_df, 2, function(x){
      naRow_v <- which(is.na(x))
      med_v <- median(x, na.rm = T)
      x[naRow_v] <- med_v
      return(x)
    })
  } # fi
  
  ## Add and subtract expressions (unless doing single-gene)
  if (ncol(sub_df) > 1){
    calc_mat <- addSubCols(sub_df)
  } else {
    calc_mat <- sub_df
  }
  comboCols_v <- colnames(calc_mat)
  
  ## Skip too few measurements
  naRow_v <- which(apply(calc_mat, 1, function(x) length(which(is.na(x)))) == ncol(calc_mat))
  if (length(naRow_v) > nrow(sub_df)*.8) {
    print(sprintf("Skipping combo %s. Fewer than 20 pct of patients had full measurements.", paste(genes_v, collapse = "_")))
    return(rep(NA, ncol(calc_mat)))
  }
  
  ## Divide patients into groups
  if (zScore_v){
    simple_mat <- apply(calc_mat, 2, function(x){
      low_v <- which(x <= -1.96); high_v <- which(x >= 1.96)
      no_v <- which(x > -1.96 & x < 1.96)
      x[low_v] <- "sigLow"; x[high_v] <- "sigHigh"; x[no_v] <- "noChange"
      return(x)
    })
  } else {
    simple_mat <- apply(calc_mat, 2, function(x){
      med_v <- median(x, na.rm = T)
      low_v <- which(x<med_v); high_v <- which(x>=med_v)
      x[low_v] <- "low"; x[high_v] <- "hi"
      #x <- factor(x, levels = c("low", "hi"))
      return(x)
    })
  }
  
  ## Combine with clin
  final_mat <- cbind(clin_df, simple_mat)
  
  ## Create survdiff objects (one for each combination)
  sd_ls <- sapply(comboCols_v, function(x) survdiff(surv_obj ~ final_mat[,x]), simplify = F)
  
  ## Extract p-value for each survival analysis
  sdp_v <- sapply(sd_ls, function(x) pchisq(x$chisq, df=1, lower=F))
  
  ## Return survival data and p-value
  return(list("values" = final_mat, "sig" = sdp_v, "check_calc" = calc_mat))
}  

addSubCols <- function(data_df){
  #' Add/Subtract column values in a specific manner
  #' @description Given a data.frame with 2 or 3 columns, create new vectors of length=nrow(data.frame) adding and subtracting
  #' the columns of the data.frame in the following manner: For 3 columns: 1+2+3, 1+2-3, 1-2+3, -1+2+3 (column indices).
  #' For 2 columns: 1+2, 1-2, -1+2.
  #' @param data_df 2- or 3-column data.frame. rows = patients/samples, columns = genes to add/subtract
  #' @value 3- or 4-column matrix consisting of results of addition/subtraction of input columns.
  #' @export
  
  ## Get order of changes
  neg_v <- c(0, ncol(data_df):1)
  
  ## Calculate additions/subtractions
  out_mat <- sapply(neg_v, function(x){
    if (x == 0) {
      out_v <- rowSums(data_df)
    } else {
      curr_df <- cbind(data_df[,-x], (data_df[,x]*-1))
      out_v <- rowSums(curr_df)
    }
    return(out_v)})
  
  ## Get gene names
  cols_v <- colnames(data_df)
  
  ## Different names depending on 2 or 3 columns
  if (ncol(data_df) == 2){
    newCols_v <- c(paste(cols_v, collapse = "+"), paste(cols_v, collapse = "-"), paste(paste0("-", cols_v[1]), cols_v[2], sep = "+"))
  } else if (ncol(data_df) == 3){
    newCols_v <- c(paste(cols_v, collapse = "+"), paste(paste(cols_v[1:2], collapse="+"), cols_v[3], sep = "-"),
                   paste(cols_v[1], paste(cols_v[2:3], collapse = "+"), sep = "-"), paste0("-", paste(cols_v, collapse = "+")))
  } else {
    stop("incorrect number of columns.")
  } # fi
  
  ## Rename
  colnames(out_mat) <- newCols_v
  
  ## Return
  return(out_mat)
}



# ## Find any with p < pthresh
# torun_v <- sdp_v[which(sdp_v <= pthresh_v)]
# ## Plot, if any
# if (length(torun_v) > 0){
#   for (col_v in names(torun_v)){
#     ## Surv fit
#     curr_s <- survfit(surv_obj ~ final_mat[,col_v])
#     ## Plot
#     if (plot_v){
#       plot(curr_s, col=c("blue", "red"), frame = F, lwd = 2, main = col_v, mark.time = T,
#            xlab = "Days to Disease Progression", ylab = "Prop. Survival",
#            cex.main = 2, cex.lab = 1.5, cex.axis = 1.2)
#       ## Add median survival line
#       medians_v <- summary(curr_s)$table[,'median']
#       lines(c(0,max(medians_v)),c(0.5,0.5), lty = "dashed")
#       mapply(function(x,y) lines(rep(x,2),c(0,0.5),col=y, lty="dashed"), medians_v, c("blue","red"))
#       ## Add legends
#       max_x <- max(final_mat$Days, na.rm = T)*.95
#       legend(max_x, 1, legend = paste0("p.val = ", round(torun_v[col_v], digits = 3)), bty = "n", xjust=1, cex = 1.2)
#       legend(max_x, 0.94, legend = c(paste0("High = ", medians_v[1]), paste0("Low = ", medians_v[2])), bty = "n", xjust = 1, col=c("blue", "red"), lwd=2,  cex = 1.2)
#     } # fi
#   } # for col_v
# } # fi