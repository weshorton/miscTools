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

suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(data.table))
library(svglite)
source("~/stable_repos_11_17/WesPersonal/utilityFxns.R")
suppressMessages(library(optparse))

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
    c("-v", "--division"),
    type = "character",
    default = "quartile",
    help = "Character vector indicating how to divide patients. quartile (default): 'low' patients are below 1st quartile, 'up' patients are above 3rd.
    median: 'low' patients are below median, 'up' patients are above median. third: 'low' patients are below bottom third, 'up' patients above top third."
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
  ),
  make_option(
    c("-y", "--type"),
    type = "character",
    default = "pdf",
    help = "'svg' - create svg plot; 'pdf' - create pdf plot."
  )
)

### Parse Commandline
p <- OptionParser(usage = "%prog -i inputFile -c countCols -s sampleCol -t timeCol -d deathCol -a altMeasure -v division -o outDir -p plot -y type",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

inputFile_v <- args$inputFile
countCols_v <- args$countCols
sampleCol_v <- args$sampleCol
timeCol_v <- args$timeCol
deathCol_v <- args$deathCol
altMeasure_v <- args$altMeasure
div_v <- args$division
outDir_v <- args$outDir
plot_v <- args$plot
type_v <- args$type

### For testing
# inputFile_v <- "/Users/hortowe/projs/Sushil/2018_02_26/Xena/data/counts_and_clinical/tcga_BLCA.tsv"
# countCols_v <- "3,4,5,6,7,8"
# sampleCol_v <- 1
# timeCol_v <- 12
# deathCol_v <- 11
# altMeasure_v <- 13
# # div_v <- "third"
# # div_v <- "median"
# div_v <- "quartile"
# outDir_v <- "/Users/hortowe/projs/Sushil/2018_02_26/Xena/survival/"
# plot_v <- F
# type_v <- "pdf

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

### Make outout
outDir_v <- mkdir(outDir_v, inputName_v)

##########################
### WRANGLE INPUT DATA ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##########################
print("Begin Wrangle")
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

###
###
### Intermediate Summary: At this point, data_dt is a data.table with 1 row per patient. Columns contain sample ID, all of the genes we decided to analyze, as well
###                       as various clinical data points (most importantly days to death/follow up and event masking)
###
###

### For each gene, we need to determine if it is above or below the cut-off in order to place it in the appropriate group.
for (i in 1:length(countCols_v)){
  ## Get column
  currCol_v <- countCols_v[i]
  newCol_v <- paste(currCol_v, "grp", sep = "_")
  
  ## Get summary. 
  if (div_v %in% c("third", "Third")) {
    currSummary_v <- thirds(data_dt[[currCol_v]], na.rm = T)
  } else {
    currSummary_v <- fivenum(data_dt[[currCol_v]], na.rm = T)
  } # fi
  
  ## Get cut-offs
  if (div_v %in% c("median", "Median")) {
    low_v <- up_v <- currSummary_v[3]
  } else {
    low_v <- currSummary_v[2]
    up_v <- currSummary_v[4]
  } # fi

  ## Divide into low and high groups
  if (div_v %in% c("median", "Median")) {
    currLowIndex_v <- which(data_dt[[currCol_v]] <= low_v)
    currUpIndex_v <- which(data_dt[[currCol_v]] > up_v)
  } else {
    currLowIndex_v <- which(data_dt[[currCol_v]] <= low_v)
    currUpIndex_v <- which(data_dt[[currCol_v]] >= up_v)
  } # fi
  
  ## Create new column and update
  data_dt[[newCol_v]] <- "mid"
  data_dt[currLowIndex_v, eval(newCol_v)] <- "low"
  data_dt[currUpIndex_v, eval(newCol_v)] <- "up"
}
print("End Wrangle")

######################
### SURVIVAL MODEL ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################

### Surv arguments are generally:
### time = end time
### event = indicates whether or not an event occurred (0 or 1)
### type = what kind of masking is performed, usually not used

### Since different patients may be used for different comparisons, we have to create a unique model for each comparison
groups_mat <- combn(countCols_v, 2)
groupNames_v <- apply(groups_mat, 2, function(x) paste(x[1], x[2], sep = "_"))

### Two-gene comparisons create 4 groups (lo_lo, lo_up, up_lo, up_up)
### Want to do everything 2 times - all 4 together and lo_lo vs. up_up

expGrps_lsv <- list("all" = c("low_low", "low_up", "up_low", "up_up"),
                    "lo_up" = c("low_low", "up_up"))


finalList_lslslsv <- list()

for (j in 1:length(expGrps_lsv)){
  ## Get groups and divisions
  currExpGrpName_v <- names(expGrps_lsv)[j]
  currExpGrp_v <- expGrps_lsv[[currExpGrpName_v]]

  ## Create empty lists to hold everything
  survObj_ls <- survDiff_ls <- pVals_ls <- survFit_ls <- list()

  for (i in 1:length(groupNames_v)){
    ## Get columns and name
    col1_v <- groups_mat[1,i]; col2_v <- groups_mat[2,i]
    name_v <- groupNames_v[i]
    
    ## Get group columns
    grp1_v <- paste0(col1_v, "_grp"); grp2_v <- paste0(col2_v, "_grp")
    
    ## Make final column
    finalCol_v <- paste0(name_v, "_grp")
    data_dt[[finalCol_v]] <- paste(data_dt[[grp1_v]], data_dt[[grp2_v]], sep = "_")
    
    ## Subset
    currData_dt <- data_dt[get(finalCol_v) %in% currExpGrp_v,]

    ## Survival Object
    survObj_ls[[name_v]] <- Surv(time = currData_dt[[timeCol_v]], event = currData_dt$Event)
    
    ## Model
    survDiff_ls[[name_v]] <- survdiff(survObj_ls[[name_v]] ~ currData_dt[[finalCol_v]])
    
    ## P-value
    pVals_ls[[name_v]] <- pchisq(survDiff_ls[[name_v]]$chisq, df = 1, lower = F)
    
    ## Fit
    survFit_ls[[name_v]] <- survfit(survObj_ls[[name_v]] ~ currData_dt[[finalCol_v]])
    
  } # for i
  
  ## Add to final list
  finalList_lslslsv[[currExpGrpName_v]] <- list("survObj" = survObj_ls, "survDiff" = survDiff_ls, "pVal" = pVals_ls, "survFit" = survFit_ls)
} # for j

print("End Model")
############
### PLOT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############

for (j in 1:length(finalList_lslslsv)) {
  ## Get name and data
  compName_v <- names(finalList_lslslsv)[j]
  currComp_lslsv <- finalList_lslslsv[[compName_v]]
  print(sprintf("Currently on: %s", compName_v))
  
  ## Make new directory
  compDir_v <- mkdir(outDir_v, compName_v)
  
  ## Make plot for each comparison
  for (i in 1:length(groupNames_v)) {
    ## Get data
    currName_v <- groupNames_v[i]
    currSplitName_v <- strsplit(currName_v, split = "_")[[1]]
    currSurvObj <- currComp_lslsv$survObj[[currName_v]]
    currSurvDiff <- currComp_lslsv$survDiff[[currName_v]]
    currP <- currComp_lslsv$pVal[[currName_v]]
    currSurvFit <- currComp_lslsv$survFit[[currName_v]]
    print(sprintf("Currently on: %s", currName_v))
    
    ## Get dimension to check
    dim_v <- dim(summary(currSurvFit)$table)
    
    ## Determine options based on comp
    if (compName_v == "all") {
      colors_v <- c("darkblue", "blue", "red", "darkred")
      ys_v <- c(0.6, 1, 0.8)
      if (dim_v[1] < 4) next
    } else {
      colors_v <- c("darkblue", "darkred")
      ys_v <- c(0.75, 1, 0.87)
      if (dim_v[1] != 2) next
    } # fi
    
    ## Get Medians and records
    medians_v <- summary(currSurvFit)$table[,'median']
    names(medians_v) <- gsub("^.*=", "", names(medians_v))
    medianLegend_v <- unlist(sapply(names(medians_v), function(x) paste0(x, " = ", medians_v[[x]]), simplify = F))
    
    records_v <- summary(currSurvFit)$table[,'records']
    names(records_v) <- gsub("^.*=", "", names(records_v))
    recordLegend_v <- unlist(sapply(names(records_v), function(x) paste0(x, " = ", records_v[[x]]), simplify = F))
    
    ## Get plot name and open file
    if (plot_v) {
      name_v <- paste("tcga", inputName_v, compName_v, currName_v, "multiGeneSurvival", sep = "_")
      if (type_v == "pdf") {
        name_v <- paste0(name_v, ".pdf")
        pdf(file = file.path(compDir_v, name_v), width = 15, height = 15)
      } else if (type_v == "svg") {
        name_v <- paste0(name_v, ".svg")
        svglite(file = file.path(compDir_v, name_v), width = 15, height = 15)
      } else {
        stop("Please choose 'pdf' or 'svg'.") 
      } # fi
    } # fi
    
    ## Generate base plot
    plot(currSurvFit, col = colors_v, frame = F, lwd = 2, cex.main = 2, cex.lab = 1.5, cex.axis = 1.2,
         mark.time = T, xlab = "Days to Death", ylab = "Prop. Survival",
         main = paste0("Survival of ", inputName_v, " Patients \nBased on ", currSplitName_v[1], " and ", currSplitName_v[2], " Expression"))
    
    ## Check if median to use
    if (length(unique(medians_v)) == 1 & is.na(unique(medians_v)[1])) {showMedians_v <- F} else {showMedians_v <- T}
    
    ## Add median lines
    if (showMedians_v){
      lines(c(0,max(medians_v, na.rm = T)),c(0.5,0.5), lty = "dashed")
      mapply(function(x,y) lines(rep(x,2),c(0,0.5),col=y, lty="dashed"), medians_v, colors_v)
    }
    
    ## Add legends
    max_x <- max(currSurvObj[,1], na.rm = T) *.95
    legend(max_x, ys_v[1], legend = paste0("p.val = ", round(currP, digits = 4)), bty = "n", xjust=1, cex=1.2)
    if (showMedians_v) legend(max_x, ys_v[2], legend = medianLegend_v, bty = "n", xjust = 1, col = colors_v, lwd = 2, cex = 1.2, title = "Median Survival")
    legend(max_x, ys_v[3], legend = recordLegend_v, bty = "n", xjust = 1, col = colors_v, lwd = 2, cex = 1.2, title = "N. Patients")
    
    if (plot_v) dev.off()
  } # for i
  
} # for j
