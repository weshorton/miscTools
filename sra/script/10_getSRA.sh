#!/bin/bash

### Download data from the Sequence Read Archive (SRA) and optionally send it through STAR alignment

#################
### Constants ###
#################

MYAPPS="$REPOS/myApps"
ESEARCH="$MYAPPS/edirect/esearch"
EFETCH="$MYAPPS/edirect/efetch"

#################
### Arguments ###
#################

EXP_ID=$1         # Should be a BioProject ID number, such as PRJNA251383
RUNSRA=$2         # logical, indicating if SRA download should be run or not
REF=$3            # path to directory for reference files
OUT1=$4           # path to desired output directory for SRA download
SUB=$5

echo $EXP_ID
echo $RUNSRA
echo $REF
echo $OUT1
echo $SUB

####################
### SRA Download ###
####################

### Get database search query:        $ESEARCH -db sra -query $EXP_ID
### Get run info from query:          $EFETCH --format runinfo
### Extract run only:                 cut -d ',' -f 1
### Get the right type:               grep SRR
### Download the files:               xargs fastq-dump --split-files

if [ $RUNSRA == T ]; then

    ## Update
    echo "Beginning SRA download."

    #    $ESEARCH -db sra -query $EXP_ID | $EFETCH --format runinfo | head -5 | cut -d ',' -f 1 | grep SRR | xargs fastq-dump --split-files -O $OUT1
    ## Get IDs
    $ESEARCH -db sra -query $EXP_ID | $EFETCH --format runinfo | cut -d ',' -f 1,30 | grep SRR > $REF/sraIDs.txt

    ## Subset for desired samples using the matching geo_ids, then cut to just get sra ids
    if [ -z "$SUB" ]; then
        echo "No subset"
    else
        echo "Subset for specific IDs"
        grep -f $REF/$SUB/er_geo_ids.txt $REF/sraIDs.txt | cut -d ',' -f 1  > $REF/$SUB/er_sra_ids.txt
    fi
fi


# ######################
# ### STAR Alignment ###
# ######################


# if [ $RUNSTAR == T ]; then

#     ## Capture SRA process
#     fullProcess=`eval condor_q | grep getSRA.sh`

#     ## Check if process exists or not (-z returns true if string is empty)
#     if [ -z "$fullProcess" ]; then

# 	## Update
# 	echo "Done with SRA download. Preparing to run STAR"

# 	## Run STAR for each file
# 	for file in `ls $OUT1`; do
	    
    
