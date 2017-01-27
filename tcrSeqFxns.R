# Sort files in a directory numerically, given a regex that substitutes all but a sample number
# inputDir_v : directory containing all and only desired files
# regex_v : regular expression that removes all portions of file name except for sample number
fileSort <- function(inputDir_v, regex_v){
  files_v <- list.files(inputDir_v)
  sortedFiles_v <- files_v[order(as.numeric((gsub(regex_v, "", files_v))))]
  return(sortedFiles_v)
}

# Read in count data using filenames from metadata. Combine counts based on unique observational varible (AA. Seq. CDR3, for example)
# inputDir_v : directory containing all and only desired count files
# inputFiles_v : vector of file names from inputDir_v
# metadata_dt : data.table of metadata containing 
# fileCol_v : colname of metadata_dt containing file names
# sampleCol_v : colname of metadata_dt containing sample numbers
# countCol_v : colname of resulting data.tables that should be included in output matrix
# obsCol_v : unique observational varaiable by which to combine columns
# Read in a data.table, sort, colapse observation variable and sum count variable.
readAndUniqData <- function(inputDir_v, metadata_dt, fileCol_v = "files", sampleCol_v = "sample",
                            countCol_v, obsCol_v, index){
  # Read in data, subset to countCol and obsCol only. Set key to be obsCol
  comboData_dt <- fread(paste0(inputDir_v, metadata_dt[index, get(fileCol_v)]))
  comboData_dt <- comboData_dt[,.(get(obsCol_v), get(countCol_v))]
  colnames(comboData_dt) <- c(obsCol_v, countCol_v)
  setkeyv(comboData_dt, obsCol_v)
  # Collapse non-unique obsCol counts
  comboData_dt <- comboData_dt[,.(countCol_v = sum(get(countCol_v))), by = .(get(obsCol_v))]
  colnames(comboData_dt) <- c(obsCol_v, metadata_dt[index,get(sampleCol_v)])
  return(comboData_dt)
} # readAndUniqData()

combineCounts <- function(inputDir_v, inputFiles_v = NULL, metadata_dt, fileCol_v = 'files', sampleCol_v = "sample",
                          countCol_v, obsCol_v){
  # get files if haven't already
  if (is.null(inputFiles_v)){
    inputFiles_v <- list.files(inputDir_v)
  } # fi
  # Read in all files and subset for desired columns
  inputData_lsdt <- list()
  for (i in 1:metadata_dt[,.N]){
    inputData_lsdt[[i]] <- readAndUniqData(inputDir_v, index = i, metadata_dt, fileCol_v, sampleCol_v,
                                           countCol_v, obsCol_v)
  } # for 
  # Outer join
  comboData_dt <- inputData_lsdt[[1]]
  for (i in 2:length(inputData_lsdt)){
    comboData_dt <- merge(x = comboData_dt, y = inputData_lsdt[[i]], by = obsCol_v, all = T)
  }
  return(comboData_dt)
} # combineCounts
