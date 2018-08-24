#!/bin/bash

### Download data from the Sequence Read Archive (SRA) and optionally send it through STAR alignment

#################
### Constants ###
#################

ESEARCH="$BIOCODERS/Applications/edirect/esearch"
EFETCH="$BIOCODERS/Applications/edirect/efetch"
FASTQDUMP="$BIOCODERS/Applications/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump"

#################
### Arguments ###
#################

FASTQ_DIR=$1            # Path to directory to hold downloaded fastq files
SRA_ID=$2               # Single ID of SRR file to download

######################
### FASTQ DOWNLOAD ###
######################

printf "Downloading ID: %s\n" $SRA_ID

cmd="$FASTQDUMP --split-files -O $FASTQ_DIR $SRA_ID"

echo $cmd
eval $cmd

