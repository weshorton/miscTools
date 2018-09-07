
######################################
### Convert BAM file to Data Table ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################

# Input: BAM file has a list of list structure. Elements in second list are columns in original BAM file
# Output: Collapsed list of list, where second-list elements are now columns in a data.table

bam.to.dt <- function(bam.file){
  # Collapse List of list to just 1 list
  my.collapse.bam <- lapply(names(bam.file[[1]]), function(elt) {do.call(c, unname(lapply(bam.file, "[[", elt)))})
  # Add names back
  names(my.collapse.bam) <- names(bam.file[[1]])
  # Turn into data.frame and then convert to data.table
  my.bam.dt <- as.data.table(do.call("DataFrame", my.collapse.bam))
  return(my.bam.dt)
} # bam.to.dt

###########################
### Update CIGAR String ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###########################

# Input: BAM file in data.table format, with original CIGAR column that has match, indel, soft clip, etc.
# Output: Same as input, except CIGAR column has been modified to be the sum of all matches in original string (everything else is removed)
# TO DO: currently only sums matches. Could add a third argument to specify a different option to sum (S, H, I, or something)

update.cigar <- function(bam.dt){
  # Extract cigar string from input data.table
  my.cigars <- bam.dt$cigar
  # Remove Insertions (I), Deletions (D), Soft-clipp (S), and Hard-clip (H)
  my.matches <- gsub("\\d+S|\\d+D|\\d+H|\\d+I", '', my.cigars)
  # Remove M, leaving just the number
  my.matches <- gsub("M", " ", my.matches)
  
  # This part is a little complicated, but don't know any other way to do it. Right now, each element of the
  # vector is a single character vector that can contain "\\d<space>" or "\\d<space>\\d<space>". Example:
  # my.matches
  # [1] "124 "     "74 "     "27 6 5 "    "102 20 53 20"
  # Four is the longest I've observed, but this solution should work regardless of length.
  # Apply over each element: split element by spaces, unlist into vector, change to numeric, take the sum.
  # output should be a vector, where elements have been "collapsed", from above example:
  # my.matches
  # [1] 124     74     38     195
  my.matches.sum <- sapply(my.matches, function(x) 
    sum(as.numeric(unlist(strsplit(x, split = ' ')))), USE.NAMES = F)
  # Update cigar string
  bam.dt$cigar <- my.matches.sum
  return(bam.dt)
} # update.cigar

############################
### Split Bam Data Table ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################

# Input: 
    # dt.list - List of BAM data.tables
    # col.to.sub - column from data.table that will be compared against value.to.sub. in order to create split
    # value.to.sub - value used for splitting data.table
    # sub.labels - become the names of the data.tables in the output list
    # comparison - type of comparison used
# Output: list of BAM data.tables
# TO DO: this is extremely complex...any way to make it better??

split.bam <- function(dt.list, col.to.sub, value.to.sub, sub.labels, comparison){
  # Create empty list for output
  subset.list <- list()
  for (i in 1:length(dt.list)){
    # Get a data table from the list
    curr.dt <- dt.list[[i]]
    # Evaluate specified column for matching (positive)
    # and not-matching (negative) of "value.to.sub"
    # change comparison accordingly
    if (comparison == "in"){
      # this one is currently simple, because proper pairs require both reads to be properly paired, so there is no
      # "split" or "between". Ideally, this portion could handle those cases, in case I wanted to subset by different flag values
      curr.positive <- curr.dt[curr.dt[,get(col.to.sub)] %in% value.to.sub]
      curr.negative <- curr.dt[!(curr.dt[,get(col.to.sub)] %in% value.to.sub)]
      curr.between <- NULL
    } else if (comparison == "glt"){
      # This is the more complicated one. I specify a range to take values greater than range, less than range, or within range.
      # CIGAR strings can vary between reads in a pair, so some may be split with one value greater than range and one value less than range.
      # Get those above and below range
      curr.positive <- curr.dt[curr.dt[,get(col.to.sub)] >= value.to.sub[2]]
      curr.negative <- curr.dt[curr.dt[,get(col.to.sub)] < value.to.sub[1]]
      # Get those within range
      curr.between <- curr.dt[curr.dt[,get(col.to.sub)] >= value.to.sub[1] & 
                                curr.dt[,get(col.to.sub)] < value.to.sub[2]]
      # extract IDs from those within range
      curr.between.ids <- curr.between$qname
      # Use IDs to take pairs that are above/below range and copy into the between group. Result is that the between
      # group holds read pairs where at least one read is within the range.
      curr.between <- rbind(curr.between, curr.positive[curr.positive$qname %in% curr.between.ids],
                            curr.negative[curr.negative$qname %in% curr.between.ids])
      # Now have to remove those reads from the above/below groups
      setkey(curr.between, `qname`)
      curr.positive <- curr.positive[!(curr.positive$qname %in% curr.between.ids)]
      curr.negative <- curr.negative[!(curr.negative$qname %in% curr.between.ids)]
    } else if (comparison == "equals"){
      # another simple one, but would be nice if this could be more complex as well. Could chane to be more like "glt"
      # where you can take a range of insert values.
      curr.positive <- curr.dt[abs(curr.dt[,get(col.to.sub)]) == value.to.sub]
      curr.negative <- curr.dt[abs(curr.dt[,get(col.to.sub)]) != value.to.sub]
      curr.between <- NULL
    } else {stop("incorrect comparison value, please choose 'in', 'glt', or 'equals'")}
    
    # Check to see if some read pairs go into different groups
    # If there is overlap between read pairs, we want to eliminate it. If a read pair has one read in positive, and one read in 
    # negative, we remove them from both and put in a new distinction - "split"
    if (length(curr.positive[curr.positive$qname %in% curr.negative$qname]$qname) != 0){
      shared.names <- curr.positive[curr.positive$qname %in% curr.negative$qname]$qname
      pure.positive <- curr.positive[!(curr.positive$qname %in% shared.names)]
      pure.negative <- curr.negative[!(curr.negative$qname %in% shared.names)]
      split.pos.neg <- curr.dt[curr.dt$qname %in% shared.names]
      
      # The rest is different depending on if there is between or not
      if (!is.null(curr.between)){
        
        # Check to see that subset worked properly. The sum of positive, negative, split, and between should
        # equal the original data.table
        if ((length(pure.positive$qname) + length(pure.negative$qname) + length(split.pos.neg$qname)) +
            length(curr.between$qname) != length(curr.dt$qname))  {
          stop("Error in splitting data using comparison ", comparison)
        } # check if correct
        
        # Paste together name from inputs
        curr.pos.label <- paste(names(dt.list)[i], sub.labels[1], "pure", sep = '_')
        curr.neg.label <- paste(names(dt.list)[i], sub.labels[2], "pure", sep = '_')
        curr.split.label <- paste(names(dt.list)[i], "split", sep = '_')
        curr.between.label <- paste(names(dt.list)[i], "between", sep = "_")
        
        # Add to list
        subset.list[[curr.between.label]] <- curr.between
        subset.list[[curr.pos.label]] <- pure.positive
        subset.list[[curr.neg.label]] <- pure.negative
        subset.list[[curr.split.label]] <- split.pos.neg
      } else { # if there is no between dt
        # Check to see that subset worked properly. Sum of positive, negative, and split should equal original dt.
        if ((length(pure.positive$qname) + length(pure.negative$qname) + length(split.pos.neg$qname)) != 
            length(curr.dt$qname))  {
          stop("Error in splitting data using comparison ", comparison)
        } # check if correct
        # If subset is correct, add to list.
        subset.list[[curr.pos.label]] <- pure.positive
        subset.list[[curr.neg.label]] <- pure.negative
        subset.list[[curr.split.label]] <- split.pos.neg
      } # If we have split and between
    } else { # This is without split, but still might have between
      if (!is.null(curr.between)){
        # positive + negative + between should equal original
        if ((length(curr.positive$qname) + length(curr.negative$qname) + length(curr.between$qname)) !=
            length(curr.dt$qname)){
          stop("Error in subsetting data using comparison ", comparison)
        } # fi
        
        # Paste together name from inputs
        curr.pos.label <- paste(names(dt.list)[i], sub.labels[1], sep = '_')
        curr.neg.label <- paste(names(dt.list)[i], sub.labels[2], sep = '_')
        curr.between.label <- paste(names(dt.list)[i], "between", sep = "_")
        
        # Add to list
        subset.list[[curr.pos.label]] <- curr.positive
        subset.list[[curr.neg.label]] <- curr.negative
        subset.list[[curr.between.label]] <- curr.between
      } else { # fi (curr.between or not)
        # Paste together name from inputs
        curr.pos.label <- paste(names(dt.list)[i], sub.labels[1], sep = '_')
        curr.neg.label <- paste(names(dt.list)[i], sub.labels[2], sep = '_')
        # Add to list
        subset.list[[curr.pos.label]] <- curr.positive
        subset.list[[curr.neg.label]] <- curr.negative
      } 
    } # fi (split output or don't, depending on results)
  } # for i
  return(subset.list)
} # split.bam.dt

# Subset based on multiple criteria
# primary.bam.dt is the dt from the bam file filtered for primary reads and sorted by qname
# raw.bam.dt is the dt from the bam file with no filter, sorted by qname
subset.bam.file <- function(primary.bam.dt, raw.bam.dt){
  ###
  ### Split bam file into those that have secondary alignments and those that do not
  ###
  # First get count of alignments to each read from raw file. If there are 3 alignments for a read pair, there is a secondary alignment
  cat("Dividing by secondary alignment or not\n")
  count.alignments <- count(raw.bam.dt$qname)
  have.2o.alignments <- count.alignments[count.alignments$freq >= 3,]
  have.2o.alignments <- as.character(have.2o.alignments$x)
  # Perform split
  no.secondary <- primary.bam.dt[!(primary.bam.dt$qname %in% have.2o.alignments)]
  yes.secondary <- primary.bam.dt[(primary.bam.dt$qname %in% have.2o.alignments)]
  # Add subsetted dt's to a list
  secondary.list <- list("no" = no.secondary, "yes" = yes.secondary)
  
  ###
  ### Split based on proper/improper pairs
  ###
  cat("Dividing by proper vs improper pairs\n")
  proper.improper <- split.bam(secondary.list, "flag", c(99, 147, 83, 163), c("proper", "improper"), "in")
  
  ###
  ### Split based on cigar < 100 and > 130
  ###
  cat("Dividing by cigar string values\n")
  cigar.100.130 <- split.bam(proper.improper, "cigar", c(100,130), c("above", "below", "between"), "glt")
  
  ###
  ### Split based on |insert length| == 200
  ###
  # cat("Dividing by insert length size\n")
  # ins.len.200 <- split.bam(cigar.100.130, "isize", 200, c("equal", "unequal"), "equals")
  
  ###
  ### Combine together
  ###
  total.subset <- list("secondary" = secondary.list, "prop.imp" = proper.improper, "cigar" = cigar.100.130) #, 
  #                      "ins.len" = ins.len.200)
} #subset.bam.file(primary.bam.dt, raw.bam.dt)

# Function to test subsets to see which are empty and also search for spike 9mers in those that are present.
# subset.list is the output of subset.bam.file(). It is a list containing lists of data.dts of different
# subsets of the original primary bam dt.
test.empty.spike.search <- function(subset.list){
  # Multiple lists within input list. Need to check all lists
  for (i in 1:length(subset.list)){
    # Empty vector for names of empty dts
    to.remove <- NULL
    # Multiple data.tables within each sublist.
    for (j in 1:length(subset.list[[i]])){
      # Check if the table is empty and add its name to the remove list if it is
      if (length(subset.list[[i]][[j]]$qname) == 0){
        to.remove <- c(to.remove, names(subset.list[[i]])[j])
      } else {
        # If it's not empty, check for both of the 9mer spike barcodes, allowing one mismatch.
        # Add results together and place all three in new columns in the dt.
        subset.list[[i]][[j]]$first.mer <- vcountPattern("TAAGTCGAC", subset.list[[i]][[j]]$seq, 
                                                         max.mismatch = 1)
        subset.list[[i]][[j]]$second.mer <- vcountPattern("GTCGACTTA", subset.list[[i]][[j]]$seq, 
                                                          max.mismatch = 1)
        subset.list[[i]][[j]]$total.mer <-
          subset.list[[i]][[j]]$first.mer + subset.list[[i]][[j]]$second.mer
      } # fi
    } # for j
    # Print all of the names of the empty data.tables in that sublist
    print(to.remove)
    # Get Ids
    #subset.ids <- names(subset.list[[i]]) %in% to.remove
    # remove them from the list.
    #subset.list[[i]][subset.ids] <- NULL
  } # for i 
  return(subset.list)
} # test.empty.spike.search