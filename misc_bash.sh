###
### Miscellaneous Bash Utilities
###



### Go through a fastq file, extract all header lines and add specific string, based on header type.
### Used with synthetic sequences in TCRseq project
### Tags: if, else, awk, OFS, fastq

awk -F '\t' -v OFS='' '{if (($1 ~ "^@") && ($1 !~ "^@NS5")) print $0,"1:11101:",NR-1,":",NR+1, " 1:N:0:ACCGAAAC"; else print $0 }' <file> | head
