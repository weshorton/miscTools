# Home for general utility functions

###
### makeDir ##################################################################################################################
###

mkdir <- function(base_dir_v, 
                  new_dir_v){
    #' Creates new directory in which to write files
    #' @description
    #' Given a base directory and string, will check if specified directory exits, and make it if not.
    #' @param base_dir_v Character string. Relative or absolute path to directory that will hold new directory
    #' @param new_dir_v Character string. Name of new directory.
    #' @return Character string of path to new directory. Makes directory in file system.
    #' @examples 
    #' makeDir("~/", "makeDir_test")
    #' @export
    
  # Concatenate to final path
  temp_dir_v <- file.path(base_dir_v, new_dir_v)
  # Add trailing slash, if absent
  if (substring(temp_dir_v, nchar(temp_dir_v)) != "/") {
    temp_dir_v <- paste0(temp_dir_v, "/")
  } # fi
  # Make if doesn't already exist
  if(dir.exists(temp_dir_v)){
    return(temp_dir_v)
  } else {
    dir.create(temp_dir_v)
    # Return string of path
    return(temp_dir_v)
  } # fi
} # makeDir

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
  # Need to get the first column and all the columns with samples
  # Sometimes last column is not updated with sample1.i and is just sample1
  allCols_v <- c("V2", countCols_v)
  altAllCols_v <- c("V2", countCols_v[1:length(countCols_v)-1], "sample1")
  
  mageckMergeCounts_dt <- tryCatch(
    {outDt <- eval(as.name(name_v))[,mget(allCols_v)]}, 
    error = function(e){outDt <- eval(as.name(name_v))[,mget(altAllCols_v)]
  })
  
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


###
### Rd to markdown #########################################################################################################
###

# from: http://stackoverflow.com/questions/27451937/render-package-documentation-into-github-wiki

rd_to_markdown <- function(rd) {
  html <-  paste(rd, ".html", sep = "")
  tools::Rd2HTML(tools::parse_Rd(rd), out = html)
  system( paste("pandoc -s -r html ", html, " -o ", rd, ".text", sep=""))
  unlink(html)
}


###
### returnSessionInfo #########################################################################################################
###

returnSessionInfo <- function(out_dir_v = NULL, load_dirs_v = NULL, args_lsv = NULL){
    #' Output session information
    #' @description Print general R session info, current date/time, and loaded git repo versions to stdout or file. Note
    #' this uses sessionInfo() from base R, but an argument could be made to use session_info() from devtools.
    #' @param out_dir_v optional directory to write information to. Default is to print to stdout
    #' @param load_dirs_v vector of directory paths that have been sourced in script. Note, if using the utility function
    #' sourceDir, these will be the same paths.
    #' @param args_lsv optional list of arguments used as input. Could be from optparse (easiest), or can hand-make a list of arguments.
    #' @value multiple lines of text printed to stdout or written to file
    #' @export

    ## Close all connections
    closeAllConnections()
    
    ## Get information
    info_v <- sessionInfo()
    date_v <- date()
    
    if (!is.null(load_dirs_v)){
        
        hashes_v <- sapply(load_dirs_v, function(x){
            setwd(x)
            system("git rev-parse --short HEAD", intern = T)
        })
        
        git_hashes_df <- as.data.frame(hashes_v)
        
    } else {
        
        git_hashes_df <- "No sourced repos"
        
    } # fi

    if (is.null(args_lsv)) args_lsv <- "No arguments given."
        

    ## Default: print to stdout
    if (is.null(out_dir_v)){
        
        cat("Session Info:\n\n"); print(info_v)
        cat("\nDate of Run:\n"); print(date_v)
        cat("\nGit repo commits:\n"); print(git_hashes_df)
        cat("\nArguments Passed:\n"); print(args_lsv)
        
    } else {
        
        ## Option: write to file
        outData_v <- paste(strsplit(as.character(Sys.time()), split = ' ')[[1]][1:2], collapse = "_")
        file_v <- file.path(out_dir_v, paste0(outData_v, "_session_info.txt"))
        file_conn <- file(file_v)
        writeLines(c("Session Info:\n", capture.output(info_v)), file_conn)
        file_conn <- file(file_v, "a")
        write(c("\nDate of Run:\n", date_v), file_conn, append = T)
        file_conn <- file(file_v, "a")
        write(c("\nGit repo commits:\n", capture.output(git_hashes_df)), file_conn, append = T)
        file_conn <- file(file_v, "a")
        write(c("\nArguments Passed:\n", capture.output(args_lsv)), file_conn, append = T)
        close(file_conn)
    } # fi

    ## Close connections
    closeAllConnections()
    
} # returnSessionInfo


###
### splitChar #################################################################################################################
###

splitChar <- function(vector_v){
    #' Split character string
    #' @description Split character string into vector where length(output) == nchar(input)
    #' @param vector_v any character vector that needs to be separated into its component parts
    #' @value another vector of length == nchar(input)
    #' @export

    out_vector_v <- unlist(strsplit(vector_v, split = ''), use.names = F)

    return(out_vector_v)
} # splitChar


###
### Detach Packages ##########################################################################################################
###

### Taken from Stack Overflow - https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r

detachAllPackages <- function() {
    #' Clear loaded packages
    #' @description Clear all loaded packages (excluding base). Useful for determining if a script has the appropriate packages loaded
    #' @value None
    #' @export

    ## Define basic packages to ignore
    basePackages_v <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")

    ## Extract attached packages from package/object list
    allPackages_v <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]

    ## Remove base packages
    removePackages_v <- setdiff(allPackages_v, basePackages_v)

    if (length(removePackages_v) > 0){
        for (package in removePackages_v) {
            detach(package, character.only=TRUE)
        } # for
    } # fi
} # detachAllPackages

###
### Merge Multiple Data.Tables ##########################################################################################################
###

mergeDTs <- function(data_lsdt, mergeCol_v, keepCol_v = NULL, ...) {
    #' Merge many data.tables together
    #' @description Take many data.tables and merge on and ID column, extracting a single column from each data.table as the column of interest
    #' @param data_lsdt list of data.tables to merge
    #' @param mergeCol_v which column from all of the data.tables to use to merge
    #' @param keepCol_v which column from all of the data.tables to use as the column of interest. If NULL, use all columns
    #' @param ... extra parameters passed to merge
    #' @value data.table with ncol == length(data_lsdt) + 1. Column names are names of list, or defaults to V1, V2,...
    #' @export

    ## Grab extra arguments
    extraParams_lsv <- list(...)

    ## Handle extra arguments
    if (!is.null(extraParams_lsv$all)){
        all_v <- extraParams_lsv$all
    } else {
        all_v <- T
    } # fi

    if (!is.null(extraParams_lsv$sort)){
        sort_v <- extraParams_lsv$sort
    } else {
        sort_v <- F
    } # fi

    ## If keepCol_v is NULL, grab all other columns
    if (is.null(keepCol_v)){
        keepCol_v <- colnames(data_lsdt[[1]])[-which(colnames(data_lsdt[[1]]) %in% mergeCol_v)]
    } # fi

    ## Create initial table by extracting the 2 columns of interest from the rest
    merge_dt <- data_lsdt[[1]][,mget(c(mergeCol_v, keepCol_v))]

    ## Create initial column names (first check if list has names and add if not)
    if (is.null(names(data_lsdt))) {
        names_v <- paste("V", 1:length(data_lsdt))
        names(data_lsdt) <- names_v
    } # fi

    if (length(keepCol_v) > 1){
        colNames_v <- c(mergeCol_v, paste(names(data_lsdt)[1], keepCol_v, sep = "_"))
    } else {
        colNames_v <- c(mergeCol_v, names(data_lsdt)[1])
    } # fi

    for (i in 2:length(data_lsdt)) {
        merge_dt <- merge(merge_dt,
                          data_lsdt[[i]][,mget(c(mergeCol_v, keepCol_v))],
                          by = mergeCol_v,
                          all = all_v, sort = sort_v)
        ## Update column names
        if (length(keepCol_v) > 1){
            colNames_v <- c(colNames_v, paste(names(data_lsdt)[i], keepCol_v, sep = "_"))
        } else {
            colNames_v <- c(colNames_v, names(data_lsdt)[i])
        } # fi
        
        ## Rename columns
        colnames(merge_dt) <- colNames_v
    } # for
    return(merge_dt)
} # mergeDTs


###
### Geometric Mean ##########################################################################################################
###

geomMean <- function(counts_v) {
  #' Find geometric mean of vector of values
  #' @param counts_v numeric vector of values
  #' @value geometric mean of vector
  #' @export
   
  ## Get product of vector
  prod_v <- prod(counts_v)
  
  ## Get length of vector
  length_v <- length(counts_v)
  
  ## Apply formula; geomMean = (X1*X2*X3*...Xn)^(1/n)
  geomMean_v <- prod_v^(1/length_v)
  return(geomMean_v)
} # geomMean


###
### 45 degree column names for pheatmap (taken from https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps)
###

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

### Have to add the following in the script that's using pheatmap in order for it to work:
# assignInNamespace(x="draw_colnames", value="draw_colnames_45",
#                  ns=asNamespace("pheatmap"))

###
### NA to O Data.table #########################################################################################################
###

naTo0 <- function(data_dt, cols_v){
  if (is.numeric(cols_v)){
    ## By column number
    for (j in seq_len(ncol(data_dt))) set(data_dt,which(is.na(data_dt[[j]])),j,0)
  } else {
    ## By column names
    for (j in names(data_dt)) set(data_dt, which(is.na(data_dt[[j]])),j,0)
  } # fi
  return(data_dt)
} # naTo0


###
### Standard ggplot theme ###############
###

my_theme <- theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
