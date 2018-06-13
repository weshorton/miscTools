#!/usr/bin/Rscript

###########################################################################################################
################################## SIMPLE KAPLAN MEIER SURVIVAL ESTIMATE ##################################
###########################################################################################################

### PURPOSE - Given expression data and clinical information, generate a simple KM survival curve and p-value

### INPUT - combined gene expression and clinical data

### OUTPUT - KM survival curve with associated significance values

####################
### DEPENDENCIES ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

library(survival)
library(survminer)
library(data.table)
library(RColorBrewer)
library(xlsx)
source("~/stable_repos_11_17/WesPersonal/utilityFxns.R")
library(optparse)

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputFile"),
    type = "character",
    help = "File containing expression data. 
    Patients/Samples are rows and genes/analytes are columns.
    Survival/clinical data are also recorded in columns"
  ),
  make_option(
    c("-c", "--countCols"),
    type = "character",
    help = "Range of numbers that are the count columns of the genes/analytes. Comma-sep, no space."
  ),
  make_option(
    c("-s", "--sampleCol"),
    type = "numeric",
    default = 1,
    help = "Column index of sample name. Defaults to 1st column. Usually called 'sample'."
  ),
  make_option(
    c("-t", "--timeCol"),
    type = "numeric",
    help = "Column index of column that will be used for survival time. Either days to death, time to recurrence, or something similar. Usually called 'days_to_death'."
  ),
  make_option(
    c('-d', "--deathCol"),
    type = "numeric",
    help = "Column index of column that indicates death status. Values should be 'alive' and 'dead'. Usually called 'vital_status'."
  ),
  make_option(
    c('-a', "--altMeasure"),
    type = "numeric",
    help = "Column index of alternative measurement. Used for alive patients that don't have a date of death. Usually called 'days_to_last_follow_up'."
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "directory to write output."
  ),
  make_option(
    c("-p", "--plot"),
    type = "logical",
    default = T,
    help = "T - write plots to output directory. F - no output."
  )
)

### Parse Commandline
p <- OptionParser(usage = "%prog -i inputFile -c countCols -s sampleCol -t timeCol -d deathCol -a altMeasure -o outDir -p plot",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

inputFile_v <- args$inputFile
countCols_v <- args$countCols
sampleCol_v <- args$sampleCol
timeCol_v <- args$timeCol
deathCol_v <- args$deathCol
altMeasure_v <- args$altMeasure
outDir_v <- args$outDir
plot_v <- args$plot

### For testing
# inputFile_v <- "/Users/hortowe/projs/Sushil/2018_02_26/Xena/data/counts_and_clinical/tcga_BLCA.tsv"
# countCols_v <- "3,4,5,6,7,8"
# sampleCol_v <- 1
# timeCol_v <- 12
# deathCol_v <- 11
# altMeasure_v <- 13
# outDir_v <- "/Users/hortowe/projs/Sushil/2018_02_26/Xena/survival/"
# plot_v <- T

###################
### PRE-PROCESS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Handle multi-args
countCols_v <- as.numeric(splitComma(countCols_v))

### Get data
data_dt <- fread(inputFile_v)

### Get columns
countCols_v <- colnames(data_dt)[countCols_v]
sampleCol_v <- colnames(data_dt)[sampleCol_v]
timeCol_v <- colnames(data_dt)[timeCol_v]
deathCol_v <- colnames(data_dt)[deathCol_v]
altMeasure_v <- colnames(data_dt)[altMeasure_v]

### Get cancer name
inputName_v <- basename(inputFile_v)
inputName_v <- gsub("tcga_|.tsv", "", inputName_v)

##########################
### WRANGLE INPUT DATA ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########################

### Turn timeCol into numeric and turn non-measurements into NA
### So far non-measurements have been displayed as "--", will update this as that changes
data_dt[get(timeCol_v) == "--", eval(timeCol_v) := NA]
data_dt[get(altMeasure_v) == "--", eval(altMeasure_v) := NA]
data_dt[[timeCol_v]] <- as.numeric(data_dt[[timeCol_v]])
data_dt[[altMeasure_v]] <- as.numeric(data_dt[[altMeasure_v]])

### Add a column for masking
### 0 is mask and 1 is measure
data_dt$Event <- 1
data_dt[get(deathCol_v) == "alive", "Event" := 0]

### If no Death time, have to replace with days to last follow up
data_dt[is.na(get(timeCol_v)), eval(timeCol_v) := get(altMeasure_v)]

### Split into list of individual data.tables
groupData_lsdt <- sapply(countCols_v, function(x) {
  y <- data_dt[,mget(c(sampleCol_v, x, timeCol_v, "Event"))]
  return(y)
}, simplify = F)

for (i in 1:length(countCols_v)){
  ## Get column
  currCol_v <- countCols_v[i]
  currData_dt <- groupData_lsdt[[currCol_v]]
  
  ## Get summary and extract quartiles
  currSummary_v <- summary(currData_dt[[currCol_v]])
  low_v <- currSummary_v[2]; up_v <- currSummary_v[5]
  
  ## Divide expression values into group
  lowIndex_v <- which(currData_dt[[currCol_v]] < low_v)
  #midIndex_v <- which(groupData_dt[[currCol_v]] >= low_v & groupData_dt[[currCol_v]] < up_v)
  upIndex_v <- which(currData_dt[[currCol_v]] > up_v)
  
  ## Turn column into character
  currData_dt[[currCol_v]] <- as.character(currData_dt[[currCol_v]])
  
  ## Update values
  currData_dt[lowIndex_v,eval(currCol_v)] <- "low"
  currData_dt[upIndex_v,eval(currCol_v)] <- "up"
  
  ## Subset
  currSubData_dt <- currData_dt[c(upIndex_v,lowIndex_v)]
  
  ## Add back
  groupData_lsdt[[currCol_v]] <- currSubData_dt
}

######################
### SURVIVAL MODEL ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################

### Surv arguments are generally:
### time = end time
### event = indicates whether or not an event occurred (0 or 1)
### type = what kind of masking is performed, usually not used

survObj_ls <- survDiff_ls <- pVals_ls <- survFit_ls <- list()

for (i in 1:length(countCols_v)) {
  ## Get column and data
  currCol_v <- countCols_v[i]
  currData_dt <- groupData_lsdt[[currCol_v]]
  
  ## Base survival object
  survObj_ls[[currCol_v]] <- Surv(time = currData_dt[[timeCol_v]], event = currData_dt$Event)
  
  ## Model
  survDiff_ls[[currCol_v]] <- survdiff(survObj_ls[[currCol_v]] ~ currData_dt[[currCol_v]])
  
  ## P-value
  pVals_ls[[currCol_v]] <- pchisq(survDiff_ls[[currCol_v]]$chisq, df = 1, lower = F)
  
  ## Fit
  survFit_ls[[currCol_v]] <- survfit(survObj_ls[[currCol_v]] ~ currData_dt[[currCol_v]])
}

############
### PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

if (plot_v) {
  name_v <- paste("tcga", inputName_v, "survival.pdf", sep = "_")
  pdf(file = file.path(outDir_v, name_v))
}

for (i in 1:length(countCols_v)){
  ## Get data
  currCol_v <- countCols_v[i]
  currData_dt <- groupData_lsdt[[currCol_v]]
  currObj <- survObj_ls[[currCol_v]]
  currDiff <- survDiff_ls[[currCol_v]]
  currP <- pVals_ls[[currCol_v]]
  currFit <- survFit_ls[[currCol_v]]
  
  ## Plot
  plot(currFit, col=c("blue", "red"), frame = F, lwd = 2, 
       main = paste0("Survival of ", inputName_v, " Patients \nBased on ", currCol_v, " Expression"),
       mark.time = T, cex.main = 2, cex.lab = 1.5, cex.axis = 1.2,
       xlab = "Days to Death", ylab = "Prop. Survival")
  
  ## Add median survival lines
  medians_v <- summary(currFit)$table[,'median']
  names(medians_v) <- gsub("^.*=", "", names(medians_v))
  lines(c(0,max(medians_v)),c(0.5,0.5), lty = "dashed")
  mapply(function(x,y) lines(rep(x,2),c(0,0.5),col=y, lty="dashed"), medians_v, c("blue", "red"))
  
  ## Get number of records
  records_v <- summary(currFit)$table[,'records']
  names(records_v) <- gsub("^.*=", "", names(records_v))
  
  ## Add legends
  max_x <- max(currData_dt[[timeCol_v]], na.rm = T) * 0.95
  legend(max_x, .5, legend = paste0("p.val = ", round(currP, digits = 4)), bty = "n", xjust=1, cex = 1.2)
  legend(max_x, 1, legend = c(paste0("Low = ", medians_v[["low"]]), paste0("High = ", medians_v[["up"]])), 
         bty = "n", xjust = 1, col=c("blue", "red"), lwd=2,  cex = 1.2, title = "Median Survival")
  legend(max_x, .73, legend = c(paste0("Low = ", records_v[["low"]]), paste0("High = ", records_v[["up"]])),
         bty  = "n", xjust = 1, col=c("blue", "red"), lwd=2, cex=1.2, title = "N. Patients")
}

if (plot_v) dev.off()