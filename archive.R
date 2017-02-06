outputData_dt <- sapply((1:metadata_dt[,.N]), function(i){
  if (i %% 10 == 0){cat(i, " ")}
  if (i == 1){
    comboData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
  } else {
    toAddData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
    # Outer join data
    comboData_dt <- merge(x = comboData_dt, y = toAddData_dt, by = obsCol_v, all = T)
  } # fi
  return(comboData_dt)})

for (i in 1:metadata_dt[,.N]){
  if (i %% 10 == 0){cat(i)}
  if (i == 1){
    comboData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
  } else {
    toAddData_dt <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v, countCol_v, obsCol_v)
    # Outer join data
    comboData_dt <- merge(x = comboData_dt, y = toAddData_dt, by = obsCol_v, all = T)
  } # fi
} # for


###
### Get mageck merge matrix
###

# Read Mageck count files and combine into matrix using sgRNA identity as "grouping" variable
getMageckMergeMatrix <- function(mageckDir_v, fileRegex_v = "^.*JT_|.count.txt"){
  require(data.table)
  # Get directory, files, and names
  mageckFiles_v <- list.files(mageckDir_v)
  mageckNames_v <- unlist(lapply(mageckFiles_v, function(x) gsub(fileRegex_v, '', x)))
  
  # Read in data
  mageckData_lsdt <- lapply(mageckFiles_v, function(x) fread(paste(mageckDir_v, x, sep = '')))
  names(mageckData_lsdt) <- mageckNames_v
  
  # Combine sgRNA and Gene into a label
  mageckLabels_lsv <-lapply(mageckData_lsdt, function(x) paste(x$sgRNA, x$Gene, sep = "_"))
  
  # Add to data.tables
  for (i in 1:length(mageckData_lsdt)){
    mageckData_lsdt[[i]] <- cbind(mageckData_lsdt[[i]], mageckLabels_lsv[[i]])
  }
  
  # Set as key
  lapply(mageckData_lsdt, function(x) setkey(x, V2))
  
  # Merge - this is not independent
  if (length(mageckData_lsdt) == 1){
    mageckMerge_dt <- mageckData_lsdt[[1]]
    mageckMergeCounts_dt <- mageckMerge_dt[,c(4,3), with = F]
  } else if (length(mageckData_lsdt) == 2){
    mageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7), with = F]
  } else if (length(mageckData_lsdt) == 3){
    tempMageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    mageckMerge_dt <- merge(tempMageckMerge_dt, mageckData_lsdt[[3]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7,10), with = F]
  } else if (length(mageckData_lsdt) == 4) {
    temp1MageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    temp2mageckMerge_dt <- merge(temp1MageckMerge_dt, mageckData_lsdt[[3]], all = T)
    mageckMerge_dt <- merge(temp2mageckMerge_dt, mageckData_lsdt[[4]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7,10,13), with = F]
  } else if (length(mageckData_lsdt) == 5) {
    temp1MageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    temp2MageckMerge_dt <- merge(temp1MageckMerge_dt, mageckData_lsdt[[3]], all = T)
    temp3MageckMerge_dt <- merge(temp2MageckMerge_dt, mageckData_lsdt[[4]], all = T)
    colnames(temp3MageckMerge_dt) <- c("V2", "sgRNA.S1", "Gene.S1", "Count.S1", "sgRNA.S2", "Gene.S2", "Count.S2", "sgRNA.S3", "Gene.S3", "Count.S3",
                                       "sgRNA.S4", "Gene.S4", "Count.S4")
    mageckMerge_dt <- merge(temp3MageckMerge_dt, mageckData_lsdt[[5]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7,10,13,16), with = F]
  } else if (length(mageckData_lsdt) == 6) {
    temp1MageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    temp2MageckMerge_dt <- merge(temp1MageckMerge_dt, mageckData_lsdt[[3]], all = T)
    temp3MageckMerge_dt <- merge(temp2MageckMerge_dt, mageckData_lsdt[[4]], all = T)
    colnames(temp3MageckMerge_dt) <- c("V2", "sgRNA.S1", "Gene.S1", "Count.S1", "sgRNA.S2", "Gene.S2", "Count.S2", "sgRNA.S3", "Gene.S3", "Count.S3",
                                       "sgRNA.S4", "Gene.S4", "Count.S4")
    temp4MageckMerge_dt <- merge(temp3MageckMerge_dt, mageckData_lsdt[[5]], all = T)
    mageckMerge_dt <- merge(temp4MageckMerge_dt, mageckData_lsdt[[6]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7,10,13,16, 19), with = F]
  } else if (length(mageckData_lsdt) == 7){
    temp1MageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    temp2MageckMerge_dt <- merge(temp1MageckMerge_dt, mageckData_lsdt[[3]], all = T)
    temp3MageckMerge_dt <- merge(temp2MageckMerge_dt, mageckData_lsdt[[4]], all = T)
    colnames(temp3MageckMerge_dt) <- c("V2", "sgRNA.S1", "Gene.S1", "Count.S1", "sgRNA.S2", "Gene.S2", "Count.S2", "sgRNA.S3", "Gene.S3", "Count.S3",
                                       "sgRNA.S4", "Gene.S4", "Count.S4")
    temp4MageckMerge_dt <- merge(temp3MageckMerge_dt, mageckData_lsdt[[5]], all = T)
    temp5MageckMerge_dt <- merge(temp4MageckMerge_dt, mageckData_lsdt[[6]], all = T)
    mageckMerge_dt <- merge(temp5MageckMerge_dt, mageckData_lsdt[[7]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7,10,13,16, 19, 22), with = F]
  } else if (length(mageckData_lsdt) == 8) {
    temp1MageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    temp2MageckMerge_dt <- merge(temp1MageckMerge_dt, mageckData_lsdt[[3]], all = T)
    temp3MageckMerge_dt <- merge(temp2MageckMerge_dt, mageckData_lsdt[[4]], all = T)
    colnames(temp3MageckMerge_dt) <- c("V2", "sgRNA.S1", "Gene.S1", "Count.S1", "sgRNA.S2", "Gene.S2", "Count.S2", "sgRNA.S3", "Gene.S3", "Count.S3",
                                       "sgRNA.S4", "Gene.S4", "Count.S4")
    temp4MageckMerge_dt <- merge(temp3MageckMerge_dt, mageckData_lsdt[[5]], all = T)
    temp5MageckMerge_dt <- merge(temp4MageckMerge_dt, mageckData_lsdt[[6]], all = T)
    temp6MageckMerge_dt <- merge(temp5MageckMerge_dt, mageckData_lsdt[[7]], all = T)
    mageckMerge_dt <- merge(temp6MageckMerge_dt, mageckData_lsdt[[8]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7,10,13,16, 19, 22, 25), with = F]
  } else if (length(mageckData_lsdt) == 12) {
    temp1MageckMerge_dt <- merge(mageckData_lsdt[[1]], mageckData_lsdt[[2]], all = T)
    temp2MageckMerge_dt <- merge(temp1MageckMerge_dt, mageckData_lsdt[[3]], all = T)
    temp3MageckMerge_dt <- merge(temp2MageckMerge_dt, mageckData_lsdt[[4]], all = T)
    colnames(temp3MageckMerge_dt) <- c("V2", "sgRNA.S1", "Gene.S1", "Count.S1", "sgRNA.S2", "Gene.S2", "Count.S2", "sgRNA.S3", "Gene.S3", "Count.S3",
                                       "sgRNA.S4", "Gene.S4", "Count.S4")
    temp4MageckMerge_dt <- merge(temp3MageckMerge_dt, mageckData_lsdt[[5]], all = T)
    temp5MageckMerge_dt <- merge(temp4MageckMerge_dt, mageckData_lsdt[[6]], all = T)
    temp6MageckMerge_dt <- merge(temp5MageckMerge_dt, mageckData_lsdt[[7]], all = T)
    colnames(temp6MageckMerge_dt) <- c("V2", "sgRNA.S1", "Gene.S1", "Count.S1", "sgRNA.S2", "Gene.S2", "Count.S2", "sgRNA.S3", "Gene.S3", "Count.S3",
                                       "sgRNA.S4", "Gene.S4", "Count.S4", "sgRNA.S5", "Gene.S5", "Count.S5", "sgRNA.S6", "Gene.S6", "Count.S6",
                                       "sgRNA.S7", "Gene.S7", "Count.S7")
    temp7MageckMerge_dt <- merge(temp6MageckMerge_dt, mageckData_lsdt[[8]], all = T)
    temp8MageckMerge_dt <- merge(temp7MageckMerge_dt, mageckData_lsdt[[9]], all = T)
    temp9MageckMerge_dt <- merge(temp8MageckMerge_dt, mageckData_lsdt[[10]], all = T)
    colnames(temp9MageckMerge_dt) <- c("V2", "sgRNA.S1", "Gene.S1", "Count.S1", "sgRNA.S2", "Gene.S2", "Count.S2", "sgRNA.S3", "Gene.S3", "Count.S3",
                                       "sgRNA.S4", "Gene.S4", "Count.S4", "sgRNA.S5", "Gene.S5", "Count.S5", "sgRNA.S6", "Gene.S6", "Count.S6",
                                       "sgRNA.S7", "Gene.S7", "Count.S7", "sgRNA.S8", "Gene.S8", "Count.S8", "sgRNA.S9", "Gene.S9", "Count.S9",
                                       "sgRNA.S10", "Gene.S10", "Count.S10")
    temp10MageckMerge_dt <- merge(temp9MageckMerge_dt, mageckData_lsdt[[11]], all = T)
    mageckMerge_dt <- merge(temp10MageckMerge_dt, mageckData_lsdt[[12]], all = T)
    mageckMergeCounts_dt <- mageckMerge_dt[,c(1,4,7,10,13,16,19,22,25, 28, 31, 34, 37)]
  }
  else {
    stop("currently not equipped to handle more than 8 samples (or 12)")
  }
  # Add column names
  colnames(mageckMergeCounts_dt) <- c("label", mageckNames_v)
  
  # Turn into matrix
  mageckMergeCounts_mat <- as.matrix(mageckMergeCounts_dt[,mget(mageckNames_v)])
  rownames(mageckMergeCounts_mat) <- mageckMergeCounts_dt$label
  return(mageckMergeCounts_mat)
} #getMageckMergeMatrix(mageckDir_v)