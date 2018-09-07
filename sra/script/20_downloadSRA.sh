#!/bin/bash

### Use fasterq-dump to download specific SRR*** IDs from a BioProject.
### This script is designed to be used after exctracting a list of IDs with 10_getSRA.sh.
### This script is submitted by the 20_sbatchDownloadSRA.submit script.

#################
### Constants ###
#################

FASTERQDUMP="$BIOCODERS/Applications/sratoolkit.2.9.2-centos_linux64/bin/fasterq-dump"

#################
### Arguments ###
#################

FASTQ_DIR=$1            # Path to directory to hold downloaded fastq files
SRA_ID=$2               # Single ID of SRR file to download

######################
### FASTQ DOWNLOAD ###
######################

printf "Downloading ID: %s\n" $SRA_ID

cmd="$FASTERQDUMP $SRA_ID -O $FASTQ_DIR --progress"
	
echo $cmd
eval $cmd

