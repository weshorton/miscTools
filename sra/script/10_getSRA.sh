#!/bin/bash

### First step in downloading data from the Sequence Read Archive (SRA)
### This tool uses the BioProject ID number for a particular dataset to grab all of the SRR*** and GSM***
### IDs from that project, to be downloaded with FASTERQDUMP later. The IDs are placed in a comma-sep file
### called sraIDs.csv (described below). The final output (and input for the next script) is final_sraIDs.csv
### and it is a single column containing only the SRR*** IDs and is optionally subset for particular IDs.

#################
### Constants ###
#################

### If you don't have access to these tools, you can download them yourself here: https://www.ncbi.nlm.nih.gov/books/NBK179288/

ESEARCH="$BIOCODERS/Applications/edirect/esearch"
EFETCH="$BIOCODERS/Applications/edirect/efetch"

#################
### Arguments ###
#################

EXP_ID=$1         # Should be a BioProject ID number, such as PRJNA251383
REF=$2            # path to directory for reference files (this is where files are output)
SUB=$3            # (optional) path to file containing the desired GSM IDs 
		  # that you would like to download (instead of the entire project)

echo $EXP_ID
echo $REF
echo $SUB

####################
### SRA Download ###
####################

### Get database search query:        $ESEARCH -db sra -query $EXP_ID
### Get run info from query:          $EFETCH --format runinfo
### Extract run only:                 cut -d ',' -f 1
### Get the right type:               grep SRR

## Update
echo "Beginning SRA download."

## Get IDs - this creates a comma-separated file (sraIDs.csv) with
	## column 1: SRR* IDs (e.g. SRR1978303)
	## column 2: GSM* IDs (e.g. GSM1660182)
	## Note that column 1 should be unique, but it is possible for column 2 to have redundancies.
	## This would indicate that a particular patient/animal has multiple fastq data files associated with it.

$ESEARCH -db sra -query $EXP_ID | $EFETCH --format runinfo | cut -d ',' -f 1,30 | grep SRR > $REF/sraIDs.csv

## Subset for desired samples using the matching geo_ids, then cut to just get sra ids

if [ -z "$SUB" ]; then
    echo "No subset"
    cut -d ',' -f 1 $REF/sraIDs.csv > $REF/final_sraIDs.csv
else
    echo "Subset for specific IDs"
    grep -f $SUB $REF/sraIDs.csv | cut -d ',' -f 1  > $REF/final_sraIDs.csv
fi


