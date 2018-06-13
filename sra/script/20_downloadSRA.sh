#!/bin/bash

### Download data from the Sequence Read Archive (SRA) and optionally send it through STAR alignment

#################
### Constants ###
#################

MYAPPS="$REPOS/myApps"
ESEARCH="$MYAPPS/edirect/esearch"
EFETCH="$MYAPPS/edirect/efetch"
FASTQDUMP="$MYAPPS/sratoolkit.2.8.2-1/bin/fastq-dump"


#################
### Arguments ###
#################

### From condor version
#var=$1
#eval `echo $var| sed -e 's/^\([^:]\{1,\}\)\:\([^:]\{1,\}\)\:\([^:]\{1,\}\)$/n=\1 mydir=\2 fil=\3/'`
#
### to do file
#myfile="$mydir/code/ref/$fil"
#mynum=`expr $n + 1`
#echo $myfile
#
#EXP_ID=`head -$mynum $myfile | tail -1`

### Slurm
mydir=$1
EXP_ID=$2

echo $EXP_ID

TEST="$FASTQDUMP --split-files -O $mydir/fastqs/ $EXP_ID"
echo $TEST
eval $TEST



# ####################
# ### SRA Download ###
# ####################

# ### Get database search query:        $ESEARCH -db sra -query $EXP_ID
# ### Get run info from query:          $EFETCH --format runinfo
# ### Extract run only:                 cut -d ',' -f 1
# ### Get the right type:               grep SRR
# ### Download the files:               xargs fastq-dump --split-files

# if [ $RUNSRA == T ]; then

#     ## Update
#     echo "Beginning SRA download."

#     #    $ESEARCH -db sra -query $EXP_ID | $EFETCH --format runinfo | head -5 | cut -d ',' -f 1 | grep SRR | xargs fastq-dump --split-files -O $OUT1
#     ## Get IDs
#     $ESEARCH -db sra -query $EXP_ID | $EFETCH --format runinfo | cut -d ',' -f 1,30 | grep SRR > $REF/sraIDs.txt

#     ## Subset for desired samples using the matching geo_ids, then cut to just get sra ids
#     grep -f $REF/er_geo_ids.txt $REF/sraIDs.txt | cut -d ',' -f 1  > $REF/er_sra_ids.txt

# fi


# # ######################
# # ### STAR Alignment ###
# # ######################


# # if [ $RUNSTAR == T ]; then

# #     ## Capture SRA process
# #     fullProcess=`eval condor_q | grep getSRA.sh`

# #     ## Check if process exists or not (-z returns true if string is empty)
# #     if [ -z "$fullProcess" ]; then

# # 	## Update
# # 	echo "Done with SRA download. Preparing to run STAR"

# # 	## Run STAR for each file
# # 	for file in `ls $OUT1`; do
	    
    
