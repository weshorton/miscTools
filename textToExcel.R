#!/usr/bin/Rscript

###
### textToExcel #################################################################################################################
###

### Convert directory of text files to Excel Workbook.

### Take directory of text files and convert to an excel workbook with one sheet per file.
### Can specify subset of files in a comma-sep list, same with preferred names. If not specified, will use
### all files in directory and file names minus extensions for the sheet names.

### Output excel file containing supplied text files, one file per sheet.

#####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(data.table)
library(optparse)
library(xlsx)
source("~/my_tool_repos/WesPersonal/utilityFxns.R")

#################
### ARGUMENTS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### Make list of options
optlist <- list(
  make_option(
    c("-d", "--dir"),
    type = "character",
    help = "Directory containing tab-sep files to combine."
  ),
  make_option(
    c("-f", "--files"),
    type = "character",
    help = "Optional parameter containing specific files inside dir_v to combine. Comma-sep."
  ),
  make_option(
    c("-w", "--wkbkName"),
    type = "character",
    help = "Name of output xlsx file."
  ),
  make_option(
    c("-s", "--sheetNames"),
    type = "character",
    help = "Optional parameter containing names to use instead of file names for sheets."
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "Optional path to output directory. If unspecified, will write to current directory."))

### Parse Command Line
p <- OptionParser(usage = "%prog -d dir -f files -w wkbkName -s sheetNames -o outDir",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

dir_v <- args$dir
files_v <- args$files
wkbkName_v <- args$wkbkName
sheetNames_v <- args$sheetNames
outDir_v <- args$outDir
  
#################
### EXECUTION ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#################

### Get files
if (is.null(files_v)) {
  files_v <- list.files(dir_v)
} else {
  files_v <- splitComma(files_v)
} # fi

### Get names
if (is.null(sheetNames_v)) {
  sheetNames_v <- tools::file_path_sans_ext(files_v)
} else {
  sheetNames_v <- splitComma(sheetNames_v)
} # fi

### Add names to files
names(files_v) <- sheetNames_v

### Prepare workbook name
if (is.null(outDir_v)) { outDir_v <- getwd() }
wkbkName_v <- file.path(path.expand(outDir_v), wkbkName_v)

### Get data
data_lsdt <- sapply(files_v, function(x) fread(file.path(dir_v, x), sep="\t", fill = T), simplify = F)

### Write files
sapply(names(data_lsdt), function(x) {
  append_v <- file.exists(wkbkName_v)
  print(x)
  write.xlsx2(data_lsdt[[x]], file = wkbkName_v, sheetName = x, row.names = F, append = append_v)
})
