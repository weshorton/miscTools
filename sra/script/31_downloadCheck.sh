#!/bin/bash

### Need to make sure all of the downloads proceded properly

### Preparation for this script
  ### After running 20_downloadSRA.sh, need to move all of the .out log files to a single directory (logDir)
  ### Need to have run 30_prepareDownloadCheck.sh to create the comparison file (refFile)

### Set ref file and log dir
refFile=$sch/publicData/code/ref/test
logDir=$sch/publicData/code/logs/sra

### Do for each output file

cd $logDir

for file in *.out; do

    ## Get write line
    wline=`tail -1 $file`

    ## Get count
    count=`echo $wline | awk -F ' ' '{print $2}'`

    ## Get name
    name=`echo $wline | awk -F ' ' '{print $5}'`

    ## Echo name
    echo $name

    ## Get ref line
    rline=`grep $name $refFile`

    ## Get ref count
    rcount=`echo $rline | awk -F ' ' '{print $4}'`

    ## Find difference
    countDiff=`expr $rcount - $count`

    ## Output, if different
    if [ $countDiff != 0 ]; then
	## Notify
	printf "Incorrect size for file %s.\nOff by: %s\n" "$name" "$countDiff"
	
	## Write out
	echo $name >> badTransfer.txt
    fi
done    
