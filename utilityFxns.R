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
###   big_label                 Similar to my_theme                      ###
###   g_legend                  extract legend from ggplot               ###
###   mhead                     print mat[1:n,1:n]                       ###
###   notice                    very easy to find print statement        ###
###   myTableGrob               custom table_grob                        ###
###   plotSpecial               plot gplot with table grobs              ###
###   mvFiles                   Move batch of files to new dir           ###
###   splitComma                Split character on the commas            ###
###   rmNARow                   Remove rows with NA in specified col     ###
###   simpleCap                 Capitalize 1st letter of each word       ###
###   anovaP                    Get p-value from ANOVA summary           ###
###   nChooseK                  Find # of unique pairwise combos         ###
###   dupCols                   Find duplicated columns (and remove)     ###
###   thirds                    Get values of min/low3rd/med/upp3rd/max  ###
###   convertDFT                Convert df to dt and vice versa          ###
###   kmPlot                    Custom Kaplain Meier plot                ###
###   quantileHeat              Heatmap with quantile color scale        ###
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
      
      ## This is new (2018-10-10) - need to make new keepCol_v if the data.tables don't have same columns
      if (!keepCol_v %in% colnames(data_lsdt[[i]])) {
        keepCol_v <- colnames(data_lsdt[[i]])[-which(colnames(data_lsdt[[i]]) %in% mergeCol_v)]
      } # fi
      
      ## Merge
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
### Angle Text theme ############################################################################################################
###

angle_x <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
angle_y <- theme(axis.text.y = element_text(angle = 45))
angle_both <- theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.text.y = element_text(angle = 45))

###
### Times New Roman #############################################################################################################
###

times <- my_theme +
  theme(plot.title = element_text(family = "Times New Roman", hjust = 0.5, size = 18),
        axis.text = element_text(family = "Times New Roman", size = 12),
        axis.title = element_text(family = "Times New Roman", size = 14),
        legend.text = element_text(family = "Times New Roman", size = 12),
        legend.title = element_text(family = "Times New Roman", size = 12),
        strip.text = element_text(family = "Times New Roman", size = 14))

###
### Big Label Times New Roman ###################################################################################################
###

bl_times <- my_theme +
  theme(plot.title = element_text(family = "Times New Roman", hjust = 0.5, size = 20),
        axis.text = element_text(family = "Times New Roman", size = 16),
        axis.title = element_text(family = "Times New Roman", size = 18),
        legend.text = element_text(family = "Times New Roman", size = 16),
        legend.title = element_text(family = "Times New Roman", size = 18),
        strip.text = element_text(family = "Times New Roman", size = 18))

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
  table_grob <- tableGrob(data_dt, rows = rep('', nrow(data_dt)), theme = ttheme_minimal())
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
   #' @export
  
   ## Logical of if row is complete or not
   completeRows_v <- complete.cases(data_dt[,..cols_v])
  
   ## Return
   if (extract_v){
     return(data_dt[completeRows_v != T,])
   } else {
     return(data_dt[completeRows_v,])
   } # fi
} # rmNARow


###
### Simple Cap ##################################################################################################################
###

simpleCap <- function(x) {
  #' Capitalize first letter of each word
  #' @description given multi-word vector, capitalize the first letter of each word.
  #' Taken from 'Examples' section of ?toupper
  #' @param x vector
  #' @value same as X, but each first letter is capitalized.
  #' @export
  
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
} # simpleCap

###
### Anova P #####################################################################################################################
###

### Extract rounded p-value from anova summary results
### THIS IS OLD! YOU SHOULD USE THE modelP FUNCTION BELOW. JUST KEEPING THIS IN CASE ANY OLD SCRIPTS USE IT.

anovaP <- function(aov, round_v = T) {
  #' Extract p-value from anova summary results
  #' @description Given the summary() results of an aov() run, extract p-value
  #' @param aov summary.aov object created by summary(aov(y ~ x))
  #' @param round_v logical. TRUE - round p-value to 3 decimals; FALSE - do not round
  #' @value numeric vector of p-value.
  #' @export
  
  aov <- aov[[1]][["Pr(>F)"]][[1]]
  aov <- ifelse(round_v, round(aov, digits = 3), aov)
  return(aov)
} # anovaP

###
### lm P #####################################################################################################################
###

lmP <- function(lm_obj, round_v = T) {
  #' Extract p-value from linear model results
  #' @description Extract the p-value estimate from a linear model output object
  #' @param lm_obj linear model object
  #' @param round_v logical. TRUE - round p-value to 3 decimals; FALSE - do not round
  #' @value numeric vector of p-value
}

###
### modelP ######################################################################################################################
###

modelP <- function(model, round_v = T) {
  #' Extract p-value from model
  #' @description Extract the p-value estimate from either a linear model (lm) or ANOVA (aov) object.
  #' @param model model object of class 'lm' for linear model, classes 'lm' and 'aov' for ANOVA.
  #' @param round_v logical. TRUE - round p-value to 3 decimals; FALSE - do not round.
  #' @value numeric vector of p-value
  #' @export
  
  ## Check class
  class_v <- class(model)
  if (length(class_v) == 2 & class_v[1] == "aov") {
    print("Supplied aov model.")
    class_v <- "aov"
  } else if (length(class_v) == 1 & class_v == "lm") {
    print("Supplied lm model.")
  } else {
    print("Supplied model is not an object of class 'aov' or 'lm'.")
  } # fi
  
  ## Take summary
  summary_v <- summary(model)
  
  ## Get p-value
  if (class_v == "aov") {
    p_v <- model[[1]][["Pr(>F)"]][[1]]
  } else {
    temp <- summary_v$fstatistic
    p_v <- pf(temp[1], temp[2], temp[3], lower.tail = F)
    attributes(p_v) <- NULL
  } # fi
  
  ## Round
  p_v <- ifelse(round_v, round(p_v, digits = 3), p_v)
  
  ## If 0, make it < 2.2e-16
  p_v <- ifelse(p_v == 0, p_v <- "< 2.2e-16", p_v)
  
  ## Output
  return(p_v)
}

###
### nChooseK ####################################################################################################################
###

nChooseK <- function(n_v, k_v) {
  #' n Choose k
  #' @description Given 'n' elements and combinations of size 'k', find number of unique combinations.
  #' @param n_v number of possible elements to create combinations
  #' @param k_v size of combinations (e.g. 3 means all unique combinations of 3 elements from n_v)
  #' @value numeric vector listing the total number of combinations
  #' @export
  
  res_v <- factorial(n_v) / (factorial(k_v) * factorial(n_v - k_v))
  return(res_v)
} # nChooseK


###
### dupCols #####################################################################################################################
###

dupCols <- function(data_dt, remove_v = T){
  #' Duplicated Columns
  #' @description Identify duplicated column names in a data.table. Print them to console and optionally remove them from data.table.
  #' Will keep the first occurrence of a column name and remove all subsequent occurrences.
  #' @param data_dt input data.table with duplicated columns
  #' @param remove_v logical, TRUE - remove columns; FALSE - do not remove, just print column names to console
  #' @value data.table same as data_dt, optionally with extra columns removed.
  #' @export
  
  ## Count occurrences of each column name
  colCounts_dt <- as.data.table(table(colnames(data_dt)))
  
  ## Extract duplicated columns
  dupCols_dt <- colCounts_dt[N > 1,]
  
  ## Get indexes to remove, for each duplicated column
  removeIndexes_v <- unlist(sapply(dupCols_dt$V1, function(x){
    y <- grep(x, colnames(data_dt))
    z <- y[2:length(y)]
    return(z)
  }, simplify = F))
  
  ## Print them
  print("Columns with more than one occurrence: ")
  print(colnames(data_dt)[removeIndexes_v])
  
  ## Remove them
  if (remove_v){
    data_dt <- data_dt[,-removeIndexes_v,with=F]
  }
  
  ## Return
  return(data_dt)
} # dupCols


###
### thirds ######################################################################################################################
###

thirds <- function(x_v, na.rm = T){
  #' Thirds-based vector division
  #' @description Instead of dividing a vector into quartiles, divide into thirds and return values at those points.
  #' @param x_v numeric vector
  #' @param na.rm TRUE (default) remove NAs. FALSE - don't remove. returns all NA
  #' @value vector containing value of min, bottom 3rd, median, top 3rd, and max.
  #' @export
  
  ## Check NAs
  isNA_v <- is.na(x_v)
  if (any(isNA_v)){
    if (na.rm) {
      x <- x[!isNA_v]
    } else {
      return(rep.int(NA,5))
    } # fi na.rm
  } # fi any()
  
  ## Sort
  x_v <- sort(x_v)
  
  ## Get length
  len_v <- length(x_v)
  
  ## Handle zero-length
  if (len_v == 0){
    rep.int(NA,5)
  } else {
    ## Get 1/4 division (len3_v*3 == len_v+1)
    len3_v <- floor((len_v+2)/1.5)/2
    ## Get divisions
    d <- c(1, len3_v, (len_v+1)/2, len_v + 1 - len3_v, len_v)
    ## Take mean of the two values surrounding division (if not an integer)
    0.5 * (x_v[floor(d)] + x_v[ceiling(d)])
  } # fi len_v == 0
} # thirds

###
### ConvertDFT ##################################################################################################################
###

convertDFT <- function(data_dft, col_v = NA, newName_v = "V1") {
  #' Convert between data.table and data.frame
  #' @description Change data.tables into data.frames with appropriate row.names or
  #' data.frames into data.tables, copying over row.names
  #' @param data_dft data in either data.table or data.frame format
  #' @param col_v if converting from dt to df, which column to use as row.names (default is 1st column)
  #' @param newName_v if converting from df to dt, what to name new column. (default is "V1")
  #' @value either a data.table or data.frame (opposite class of input)
  #' @export
  
  ## Get class
  class_v <- class(data_dft)
  
  ## Convert data.table to data.frame
  if ("data.table" %in% class_v){
    
    ## Convert
    out_dft <- as.data.frame(data_dft)
    
    ## Get column for row names
    if (is.na(col_v)) col_v <- colnames(data_dft)[1]
    
    ## Add row names
    rownames(out_dft) <- data_dft[[col_v]]
    
    ## Remove column that provided rownames
    whichCol_v <- which(colnames(data_dft) == col_v)
    out_dft <- out_dft[,-1, drop = F]
    
  ## Convert data.frame to data.table
  } else if (class_v %in% c("data.frame", "matrix")){
    
    ## Convert
    out_dft <- as.data.table(data_dft)
    
    ## Add back row.names
    out_dft[[newName_v]] <- rownames(data_dft)
    
    ## Change order
    out_dft <- out_dft[ , c(ncol(out_dft), 1:(ncol(out_dft)-1)), with = F]
  } else {
    stop("Neither 'data.table', 'data.frame', nor 'matrix' were in the class of data_dft. Please check your input data.")
  } # fi
  
  ## Return
  return(out_dft)
  
} # convertDFT
      
###
### kmPlot ######################################################################################################################
###

kmPlot <- function(fit, model, colors_v = c("blue", "red"), labels_v = c("Low" = "low", "High" = "up"), main_v = "KM Survival", 
                   xlab_v = "Overall Survival", ylab_v = "Prop. Survival", max_x, leg_y = c('med' = 1, 'rec' = 0.83, 'p' = 0.7),
                   fileName_v = NA) {
  #' Custom Kaplan Meier plot with median lines and legends
  #' @description Specific Kaplan Meier plot
  #' @param fit "survfit" object created by survfit() function
  #' @param model "survdiff" object created by survdiff() function
  #' @param colors_v vector of colors to plot. Must be same length as number of divisions in fit
  #' @param labels_v named vector of divisions in survival model. Names will be displayed on legend, values must match those of fit/model.
  #' @param main_v title for plot
  #' @param xlab_v label for X-axis. Default is "Overall Survival"
  #' @param ylab_v label for y-axis. Default is "Prop. Survival"
  #' @param max_x Max x-value from original data that created fit. Used for placing legends
  #' @param leg_y named vector of legend y-axis placements. Names are 'med', 'p', and 'rec'. Legends will only be printed if they have a leg_y value.
  #' @param fileName_v path to output file where plot will be printed. If NA, will print to stdout.
  #' @value plot to console of Kaplan Meier survival estimate
  #' @export
  
  if (!is.na(fileName_v)) pdf(file = fileName_v, width = 10, height = 10)
  ## Generate base plot
  plot(fit, col = colors_v, frame = F, lwd = 2,
       main = main_v, xlab = xlab_v, ylab = ylab_v,
       mark.time = T, cex.main = 2, cex.lab = 1.5, cex.axis = 1.2)
  
  ## Get medians and records
  medians_v <- summary(fit)$table[,'median']; names(medians_v) <- gsub("^.*=", "", names(medians_v))
  records_v <- summary(fit)$table[,'records']; names(records_v) <- gsub("^.*=", "", names(records_v))
  
  ## Add median lines
  lines(c(0, max(medians_v)), c(0.5, 0.5), lty = "dashed")
  mapply(function(x,y) lines(rep(x,2),c(0,0.5),col=y, lty="dashed"), medians_v, colors_v)
  
  ## Construct median and record legend
  if ('med' %in% names(leg_y) & 'rec' %in% names(leg_y)) {
    final_leg <- sapply(names(labels_v), function(x) {
      paste0(x, " = ", medians_v[[ labels_v[x] ]], " ; ", records_v[[ labels_v[x] ]])})
    legendTitle_v <- "Med. Surv. ; N. Records"
  } else if ('med' %in% names(leg_y) & !('rec' %in% names(leg_y))) {
    final_leg <- sapply(names(labels_v), function(x) paste(x, medians_v[[labels_v[x] ]], sep = " = "))
    legendTitle_v <- "Med. Surv."
  } else if (!('med' %in% names(leg_y)) & 'rec' %in% names(leg_y)) {
    final_leg <- sapply(names(labels_v), function(x) paste(x, records_v[[ labels_v[x] ]], sep = " = "))
    legendTitle_v <- "N. Records"
  }
  
  ## Make it
  if ('med' %in% names(leg_y) | 'rec' %in% names(leg_y)){
    legend("topright", legend = final_leg, col = colors_v, title = legendTitle_v, bty = "n", xjust = 1, cex = 1.2, lwd = 2)
  }
  # ## Add median legend
  # if ('med' %in% names(leg_y)) {
  #   med_leg <- sapply(names(labels_v), function(x) paste(x, medians_v[[labels_v[x] ]], sep = " = "))
  #   legend(max_x, leg_y['med'], legend = med_leg, col = colors_v, title = "Med. Surv.", bty = 'n', xjust = 1, cex = 1.2, lwd = 2)
  # } # fi
  # 
  # ## Add Record legend
  # if ('rec' %in% names(leg_y)) {
  #   rec_leg <- sapply(names(labels_v), function(x) paste(x, records_v[[ labels_v[x] ]], sep = " = "))
  #   legend(max_x, leg_y['rec'], legend = rec_leg, col = colors_v, title = "N. Records", bty = 'n', xjust = 1, cex = 1.2, lwd = 2)
  # } # fi
  
  ## Add pvalue legend
  if ('p' %in% names(leg_y)) {
    pval <- pchisq(model$chisq, df = 1, lower = F)
    pv_leg <- paste0("p.val = ", round(pval, digits = 4))
    legend("bottomleft", legend = pv_leg, bty = 'n', xjust = 1, cex = 1.2)
  }
  
  ## Close device
  if (!is.na(fileName_v)) graphics.off()
  
}

###
### Quantile Heat ###############################################################################################################
###

#dir_v <- dirname(sys.frame(1)$ofile)
source("~/my_tool_repos/WesPersonal/quantileHeatmap.R")

###
### Mean Pos Only ###############################################################################################################
###
  
### Answer to SO Post - https://stackoverflow.com/questions/52973875/is-there-a-simpler-way-to-calculate-mean-while-dropping-negative-values/52975841#52975841
  
meanPosOnly <- function(data, refCol_v, calcCol_v = NA, negCountName_v = "tnegcount", meanName_v = "averagev", rename_v = T) {
  #' Calculate means of positive values
  #' @description Calculate the mean value of all positive values in all rows of a data.frame, matrix, etc.
  #' @param data - data.frame, matrix, etc. Table of values
  #' @param refCol_v - character vector - Name of column(s) that will not be used in taking the mean.
  #' Some sort of reference/metadata column(s). Must be before other columns.
  #' @param calcCol_v - vector (character or numeric) -
  #' character - Name of column(s) that will be used in taking the mean. Default is NA, which will use all columns not in refCol_v.
  #' numeric - column indices of column(s) that will be used in taking the mean.
  #' @param negCountName_v - character vector - name of column that will tally number of negative values in each row
  #' @param meanName_v - character vector - name of column that will contain the resulting average of all positive values in each row
  #' @param rename_v - logical - rename the calc columns by adding ".[0-9]" where [0-9] is 1 more than currently in name
  #' @value data.frame of same dimensions as data, with 2 extra columns denoting the number of negatives in each row and the mean of all positive values.
  #' @export
  
  ## Get column indices
  if (is.na(calcCol_v[1])) {
    cols_v <- grep(paste(refCol_v, collapse = "|"), colnames(data), invert = T)
  } else if (is.character(calcCol_v)) {
    cols_v <- which(colnames(data) %in% calcCol_v)
  } else {
    cols_v <- calcCol_v
  } # fi

  ## Get numeric columns
  whichNum_v <- which(sapply(data, class) == "numeric")

  ## Get result
  out_df <- as.data.frame(t(apply(data, 1, function(x) {
    whichMean_v <- which(as.numeric(x[cols_v]) >= 0)
    num0_v <- length(cols_v) - length(whichMean_v)
    y <- mean(as.numeric(x[cols_v][whichMean_v]))
    z <- c(x, num0_v, y)
    return(z)
  })))

  ## Add names
  colnames(out_df)[c(ncol(out_df)-1,ncol(out_df))] <- c(negCountName_v, meanName_v)

  ## Fix numeric columns
  for (c_v in c(whichNum_v, ncol(out_df)-1, ncol(out_df))) out_df[,c_v] <- as.numeric(as.character(out_df[,c_v]))

  ## Fix calc names
  colNames_v <- colnames(out_df)[cols_v]
  if (rename_v) {
    colNames_v <- sapply(colNames_v, function(x) {
      y <- strsplit(x, split = "\\.")[[1]]
      z <- ifelse(is.na(y[2]),
                  paste0(y, ".1"),
                  paste0(y[1], ".", (as.numeric(y[2])+1)))
      return(z)})
    } # fi
  colnames(out_df)[cols_v] <- colNames_v

  ## Return
  return(out_df)
} # meanPosOnly