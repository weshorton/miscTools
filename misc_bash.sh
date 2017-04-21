###
### Miscellaneous Bash Utilities
###



### Go through a fastq file, extract all header lines and add specific string, based on header type.
### Used with synthetic sequences in TCRseq project
### Tags: if, else, awk, OFS, fastq

awk -F '\t' -v OFS='' '{if (($1 ~ "^@") && ($1 !~ "^@NS5")) print $0,"1:11101:",NR-1,":",NR+1, " 1:N:0:ACCGAAAC"; else print $0 }' <file> | head


### Paste two tab-sep files together
### Tags: paste, txt, tsv
paste file1 file2 > file3

### Rearrange columns
### Tags: txt, tsv, awk
awk -F '\t' -v OFS='\t' '{print $1, $3, $2}' file > newfile

### Remove trailing whitespace
### Tags: sed, whitespace, bash
sed -e 's/\s*//g' filename

### Echo statments to stderr instead of stdout
### stderr; bash; condor; stdout
echoerr() { printf "%s\n" "$*" >&2; }


### Rename multiple files based on delimiter
### bash; rename; delimiter
### If file names are: LIB170313LC_S9_R2_001.fastq.gz; and want them to be LIB170313LC_S9_S9_R2_001.fastq.gz
for i in *.gz; do mv "$i" "$(echo $i | awk -F "_" -v OFS="_" '{print $1, $2, $2, $3, $4}')"; done

