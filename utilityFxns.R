# Home for general utility functions

###
### mkdir ##################################################################################################################
###

mkdir <- function(baseDir_v, 
                  newDir_v){
  # Create new file for writing (if doesn't exist)
  # baseDir_v: some sort of file path to where new files should be written
  # newDir_v:  character vector of name of new directory to write files to
  tempDir_v <- file.path(baseDir_v, newDir_v)
  if(dir.exists(tempDir_v)){
    return(tempDir_v)
  } else {
    dir.create(tempDir_v)
    return(tempDir_v)
  } # fi
} # mkdir

###
### getMageckMergeMatrix ###################################################################################################
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
  
  # Empty variable to hold count columns
  countCols_v <- NULL
  
  # Perform merge
  for (i in 1:length(mageckData_lsdt)){
    # Get temporary name for data.table
    name_v <- paste0("temp", i)
    countCols_v <- c(countCols_v, paste0("sample1.", i))
    if (i ==1 ) {
      # Assign first data.table
      assign(name_v, mageckData_lsdt[[i]])
    } else {
      # Get first name 
      firstName_v <- paste0("temp", i-1)
      assign(name_v, merge(eval(as.name(firstName_v)), mageckData_lsdt[[i]], all = T, suffixes = c(paste0(".", i-1), paste0(".", i))))
    } # fi
  } # for i
  # Need to get the first column and all the columns with samplex
  allCols_v <- c("V2", countCols_v)
  mageckMergeCounts_dt <- eval(as.name(name_v))[,mget(allCols_v)]
  
  # Add column names
  colnames(mageckMergeCounts_dt) <- c("label", mageckNames_v)
  # Turn into matrix
  mageckMergeCounts_mat <- as.matrix(mageckMergeCounts_dt[,mget(mageckNames_v)])
  rownames(mageckMergeCounts_mat) <- mageckMergeCounts_dt$label
  return(mageckMergeCounts_mat)
} # getMageckMergeMatrix

###
### lcat ###################################################################################################################
###

# Cat information to log file rather than stdout in the console
lcat <- function(..., file_v = logFile_v){
    output_v <- c(..., "\n")
    cat(output_v, file = file_v)
    } # lcat()

###
### listToDataTable ########################################################################################################
###

# Transform a list of vectors of different lengths into a data.table
# Pad shorter vectors with NA's
listToDataTable <- function(input_lsv){
  
  # Get maximum vector length
  maxLength_v <- max(unlist(lapply(input_lsv, function(x) length(x))))
  
  # Function to extend each vector in list with NA's to match length of maximum vector
  addNAToVector <- function(listElement_v, maxLength_v){
    difference_v <- maxLength_v - length(listElement_v)
    NAToAdd_v <- rep(NA, times = difference_v)
    newElement_v <- append(listElement_v, NAToAdd_v)
    return(newElement_v)
  } # addNAToVector
  
  # Apply NA-adding function to each vector in list
  output_lsv <- lapply(input_lsv, function(x,y) addNAToVector(x, maxLength_v))
  
  # Convert to data.table
  output_dt <- as.data.table(output_lsv)
  
  # Return
  return(output_dt)
} # listToDataTable()

###
### matToDataTable #########################################################################################################
###

# Take a matrix and conver to data.table, while taking the rowNames and making them the 1st column of data.table

matToDataTable <- function(inputData_mat, colName_v){
  # Convert
  outputData_dt <- as.data.table(inputData_mat)
  # Add rownames as column
  outputData_dt[[colName_v]] <- rownames(inputData_mat)
  # Reorder columns
  outputData_dt <- outputData_dt[,c(ncol(outputData_dt), (1:(ncol(outputData_dt)-1))), with = F]
  # return
  return(outputData_dt)
} # matToDataTable

