#!/bin/sh

### Reformat sra results file downloaded from ncbi.nlm.nih.gov
### Steps to download:
  ### 1. Go to project page at GEO
  ### 2. Click on BioProject ID
  ### 3. Under "Project Data" table, Click link in 'Number of Links' column, with row value in column 'Resource Name' as 'SRA Experiments'
  ### 4. In new page, click on 'Send to' drop down (to the left of 'Manage filters' and below the search bar)
  ### 5. Choose destination -> File

### Output format:
  ### 3 Columns
    ### 1. Experiment Title - (GSM[0-9]*) 
    ### 2. SRA ID - SRR[0-9]*
    ### 2. Total Size; Mb - size in megabytes
    ### 3. Total Spots - integer
  ### Tab-separated and sorted by Experiment Title

### A few notes regarding below manipulations
  ### 1. File is comma-separated, but one of the column names also has a comma in it. 1st sed expression replaces comma with semi-colon
  ### 2. Cut to extract the 3 desired columns
  ### 3. Experiment Title column begins with the GSM ID, but has other info as well, 2nd sed expression removes extra info
  ### 4. All of the information is in double-quotes. Remove double quotes with 3rd sed expression.
  ### 5. Change from comma-sep to tab-sep
  ### 6. Remove header line
  ### 7. Sort by GSM ID

IN=$1 	 	# Resulting csv file from steps above
IN2=$2          # This is the sraIDs.csv file created by 10_getSRA.sh 
                # it needs to be sorted by GSM ID and translated from comma-sep to tab-sep
OUT=$3          # Reformatted file

### REFORMAT REFERENCE FILE
sed -E 's/"([A-Za-z ]*),( [A-Za-z]*)"/"\1;\2"/g' $IN > temp		# 1
cut -d ',' -f 2,10,12 temp > temp2					# 2
sed -E 's/:[ 0-9A-Za-z;-]*//' temp2 > temp3				# 3
sed 's/"//g' temp3 > temp4						# 4
cat temp4 | tr ',' '\t' > temp5						# 5
tail -n +2 temp5 > temp6						# 6
sort -k 1 temp6 > $OUT							# 7
rm temp*

### REFORMAT ID FILE
cat $IN2 | tr ',' '\t' > temp
sort -k 2 temp > $IN2

### JOIN FILES
join -1 1 -2 2 -o '0,2.1,1.2,1.3' $OUT $IN2 > tempJoin
mv -f tempJoin $OUT
