#!/bin/sh

### Count all mapped reads (works because 0x04 bit means read was unmapped and -F means grab everything WITHOUT that bit
samtools view -F 0x04 -c file.bam

### Count all primary alignments (0x100 bit means secondary alignment)
samtools view -F 0x100 -c file.bam

