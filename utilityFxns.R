# Home for general utility functions

library(ggplot2)
library(grid); library(gridExtra); library(gtable)

############################ Table of contents #############################
###                                                                      ###
###   mkdir                     create new directory                     ###
###   getMageckMergeMatrix      (old) merge files                        ###
###   lcat                      cat to log rather than console           ###
###   listToDataTable           xform list of vector to dt. Pad with NA  ###
###   matToDataTable            xform mat to dt. rownames(mat) = dt[,1]  ###
###   rd_to_markdown            not sure...                              ###
###   returnSessionInfo         makes detailed log of session            ###
###   splitChar                 split character string                   ###
###   detachAllPackages         clear loaded packages                    ###
###   mergeDTs                  merge dts                                ###
###   geomMean                  geometric mean of vector                 ###
###   draw_colnames_45          adjust pheatmap colnames                 ###
###   naTo0                     convert NA in a dt to 0, by column       ###
###   my_theme                  custom ggplot2 theme                     ###
###   g_legend                  extract legend from ggplot               ###
###   mhead                     print mat[1:n,1:n]                       ###
###   notice                    very easy to find print statement        ###
###   myTableGrob               custom table_grob                        ###
###   plotSpecial               plot gplot with table grobs              ###
###   mvFiles                   Move batch of files to new dir           ###
###   splitComma                Split character on the commas            ###
###   rmNARow                   Remove rows with NA in specified col     ###
###                                                                      ###
############################################################################

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
### NA to O Data.table ##########################################################################################################
###

naTo0 <- function(data_dt, cols_v){
  #' Convert NA's in a data.table to 0
  #' @description Given a set of columns, convert all NA values in those columns to 0
  #' @param data_dt data.table with NAs in some columns
  #' @param cols_v vector of either column indices or column names that have NAs to be replaced
  #' @value data.table with same dim(), but NAs are now 0 for all cols_v
  #' @export
  
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
### Standard ggplot theme #######################################################################################################
###

my_theme <- theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))

###
### Big Label  ggplot theme #####################################################################################################
###

big_label <- theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.text = element_text(size = 16),
          axis.text.y = element_text(angle = 45),
          axis.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18))


###
### Extract ggplot legend #######################################################################################################
###

### Taken from stack overflow...can't remember link

g_legend <- function(a.gplot){
  #' Extract ggplot legend
  #' @description Extract legend as separate table from ggplot object
  #' @param a.gplot a ggplot object with a legend
  #' @value a gtable object of the legend
  #' @export
  
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

###
### Special head ################################################################################################################
###

mhead <- function(data_df, n = 5){
  #' Special head for table-like objects
  #' @description Show first n rows and columns of data.table/data.frame/matrix. Good for tables with many columns
  #' @param data_df any data.table, data.frame, or matrix
  #' @param n number of rows and columns to subset by. Default is 5
  #' @value prints to console data_df[1:n,1:n]
  #' @export
    
    print(data_df[1:n,1:n])
} # mhead

###
### Notice ######################################################################################################################
###

notice <- function(statement_v) {
  #' Print easy-to-find statement
  #' @description Print a specified statement after a large ASCII "LOOK" statement.
  #' @param statement_v Some sort of character vector. Can be a paste() call, an sprintf message, etc.
  #' @value print to console
  #' @export
    
    cat("#   ### ### # #\n#   # # # # ##\n### ### ### # #\n\n")
    cat(statement_v)
    cat("\n\n\n")
} # notice


###
### Table Grob ##################################################################################################################
###

myTableGrob <- function(data_dt, title_v, fontsize_v = 14){
  #' Create custom table grob with title
  #' @description Creates table grob in format that is most common for my usage.
  #' @param data_dt Data.table that the grob will be made out of
  #' @param title_v Title for display
  #' @param fontsize_v Fontsize for title. Default is 14 (goes will with my_theme)
  #' @value gtable object
  #' @export
  
  ## Table
  table_grob <- tableGrob(data_dt, rows = rep('', nrow(data_dt)))
  ## Title
  title_grob <- textGrob(title_v, gp = gpar(fontsize = fontsize_v))
  ## Add title
  table_grob <- gtable_add_rows(table_grob, heights = grobHeight(title_grob) + unit(5,'mm'), pos = 0)
  table_grob <- gtable_add_grob(table_grob, title_grob, 1, 1, 1, ncol(table_grob))
}


###
### Plot Special ################################################################################################################
###

plotSpecial <- function(grobs_ls, loc_lsv = list("vpPlot" = c(width = 0.65, height = 1, x = 0.3, y = 0.5),
                                                 "vpLeg" = c(width = 0.25, height = 0.7, x = 0.7, y = 0.88),
                                                 "vpTab" = c(width = 0.45, height = 0.3, x = 0.7, y = 0.6)),
                        file_v = NA){
  #' Plot ggplot with extra table grob
  #' @description Plot ggplot object, legend, and at least one other information table
  #' @param grobs_ls list of grobs to plot (must be same length as loc_lsv).
  #' @param loc_lsv list of viewport location arguments. Each list element must have 4 elements in order of width, height, x, y. Must be same length as grobs_ls. vpPlot, vpLeg, and vpTab are required.
  #' @param file_v file to output plot to. Must be pdf.
  #' @details Works differently depending on if main plot is a ggplot2 object or a pheatmap (or other) object. For ggplot2 object: Minimum length is 3 objects (plot, legend, table). Table should
  #' be created using myTableGrob or similar. Legend should be made using g_legend. Defaults work just fine. For pheatmap object, the pheatmap must be saved to an object by using
  #' myObj <- grid::grid.grabExpr(pheatmap(data_dt, options)). myObj is then the first grob in grobs_ls. As of right now, I can't extract the legend from pheatmap, so grobs_ls will actually only
  #' be a length of 2 as a minimum (heatmap and table). The loc_lsv can still work with the default 3 (the vpLeg will be removed), but will also work if just two are specified. When ggplot2 is used
  #' the loc_lsv names aren't important, but when pheatmap is used, the loc_lsv names must be vpPlot and vpTab.
  #' @value print plot to viewer or as pdf
  #' @export
  
  ## Make Viewports
  vp_ls <- sapply(loc_lsv, function(z) viewport(width = z[1], height = z[2], x = z[3], y = z[4]), simplify = F)
  
  ## Subset Viewports if not ggplot
  if (!"ggplot" %in% class(grobs_ls[[1]])){
    vp_ls[["vpLeg"]] <- NULL
  }
  
  ## Open file and/or new page
  if (!is.na(file_v)) pdf(file = file_v)
  grid.newpage()
  
  ## ggplot printing
  if ("ggplot" %in% class(grobs_ls[[1]])){
    print(grobs_ls[[1]] + theme(legend.position = "none"), vp = vp_ls[[1]])
  }
  
  ## non-ggplot printing
  if (!"ggplot" %in% class(grobs_ls[[1]])){
    pushViewport(vp_ls[[1]])
    grid.draw(grobs_ls[[1]])
  }
  
  ## Plot rest
  for (i in 2:length(vp_ls)){
    upViewport(0)
    pushViewport(vp_ls[[i]])
    grid.draw(grobs_ls[[i]])
  } # for i
  if (!is.na(file_v)) dev.off()
  
}

###
### Move Logs ###################################################################################################################
###

mvFiles <- function(origDir_v, newDir_v, pattern_v = "*\\.log"){
  #' Move a batch of files to a new directory
  #' @description Originally developed to move venn diagram log output. Can be used to move any unique set of files to new directory.
  #' @param origDir_v path to original directory containing the files
  #' @param newDir_v one of two. (1) path to new directory where files should be moved to. (2) character string of new sub-directory of origDir_v to put files into.
  #' @param pattern_v some sort of regex to uniquely identify files to move.
  #' @value moves files at system level
  #' @export
  
  ## Make new directory as sub-directory, if specified
  if (length(strsplit(newDir_v, split = "/")[[1]]) == 1){
    newDir_v <- mkdir(origDir_v, newDir_v)
  } # fi
  
  ## Get files to move
  files_v <- list.files(origDir_v, pattern = pattern_v)
  
  ## Move them
  for (file_v in files_v) file.rename(from = file.path(origDir_v, file_v),
                                      to = file.path(newDir_v, file_v))
} # mvFiles

###
### Split Comma #################################################################################################################
###

splitComma <- function(string_v){
    #' Split a string along the commas in the string
    #' @description Just a concise way to split a string on the commas. Returning a vector of strings instead.
    #' @param String that will be split
    #' @value Vector of strings. One element for each section of string that was previously divided by a comma
    #' @export

    vector_v <- unlist(sapply(string_v, function(x) strsplit(x, split = ","), USE.NAMES = F))
    return(vector_v)
} # splitComma

###
### Rm NA Row ###################################################################################################################
###

rmNARow <- function(data_dt, cols_v, extract_v = F){
   #' Remove rows that have NA in certain columns
   #' @description Remove all rows that have an NA in at least one of the specified column
   #' @param data_dt data.table to have rows removed
   #' @param cols_v columns to check for NA in
   #' @param extract_v logical, TRUE - output the NA rows; FALSE (default) - output the input table with NA rows removed
   #' @value data.table either containing the offending rows (extract_v == T), or the input table with those rows removed
   #' export
  
   ## Logical of if row is complete or not
   completeRows_v <- complete.cases(data_dt[,..cols_v])
  
   ## Return
   if (extract_v){
     return(data_dt[completeRows_v != T,])
   } else {
     return(data_dt[completeRows_v,])
   } # fi
} # rmNARow