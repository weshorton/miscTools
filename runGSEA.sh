#!/bin/sh

function runGSEA ()
{
    local countData=""
    local clsData=""
    local setList=""
    local comparisons=""
    local outDir=""
    local execute=true

    local OPTIND
    local OPTARG
    local opt

    while getopts ":c:m:s:v:o:h" opt; do

	case "$opt" in
	    h)
		echo "\
		usage:
		------
		runGSEA [ -h ] [ -c countData ] [ -m ctlData ] [ -s setList ] [ -v comparisons ] [ -o outDir ]

		description:
		------------
		Run command-line version of GSEA with standard options.

		arguments:
		----------
		-h     Print this help message and exit.

		-c     Path to file containing tab-separated expression data with the following format:
		          COL1 : NAME (gene names)
		          COL2 : DESCRIPTION (can be all na, or can be some sort of description)
		          COL3-N : sample names and their expression values in the rows

		-m     Path to file containing .cls file describing phenotypes of samples. Format:
		          ROW1 : [ total \# of samples ] [ number of different ptypes ] [ 1 ]; eg. 6 2 1
		          ROW2 : [ \# ] [ TREAT1 ] [ TREAT2 ] [ ... ] [ TREATN ]; e.g. \# Veh PLX
		          ROW3 : phenotype label for each sample, 0 = TREAT1, 1 = TREAT2; e.g. 0 0 0 1 1 1

		-s     Path to file containing list of GSEA gene sets to use in analysis. One name per line.
		          e.g. h.all.v6.0.symbols.gmt

		-v     Treatment to be compared. Must be in the format of TREAT1_versus_TREAT2 from the .cls file

                -o     Path to base output directory
		"
		execute=false
		;;
	    c)
		countData="$OPTARG"
		;;
	    m)
		clsData="$OPTARG"
		;;
	    s)
		setList="$OPTARG"
		;;
	    v)
		comparison="$OPTARG"
		;;
	    o)
		outDir="$OPTARG"
		;;
	    ?)
	    echo "Error: did not recognize option, ${OPTARG}."
	    echo "Please try -h for help."
	    execute=false
	    ;;
	esac
    done

    if [[ $execute == true ]]; then
	
	## Make base output directory, if not already there
	mkdir -p $outDir

	## Run for each gene setList
	for list in `cat $setList`; do

	    echo "List: " $list
	    echo "Count: " $countData
	    echo "CLS: " $clsData#$comparison

	    ## Make final output
	    OUT=$outDir/$list/
	    mkdir -p $OUT
	    
	    java -Xmx512m -cp /Users/hortowe/jar/gsea-3.0.jar xtools.gsea.Gsea \
	    	 -res $countData \
	    	 -cls $clsData#$comparison \
     	         -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets_final/$list \
	    	 -collapse false \
	    	 -mode Max_probe \
	    	 -norm meandiv \
	    	 -nperm 1000 \
	    	 -permute phenotype \
	    	 -rnd_type no_balance \
	    	 -scoring_scheme weighted \
	    	 -rpt_label $comparison \
	         -metric Signal2Noise \
	    	 -sort real \
	    	 -order descending \
	    	 -create_gcts false \
	    	 -create_svgs false \
	    	 -include_only_symbols true \
	    	 -make_sets true \
	    	 -median false \
	    	 -num 100 \
	    	 -plot_top_x 20 \
	    	 -rnd_seed timestamp \
	    	 -save_rnd_lists false \
	    	 -set_max 500 \
	    	 -set_min 15 \
	    	 -zip_report false \
	    	 -out $OUT/ \
	    	 -gui false > $OUT/stdout_$comparison.txt
	done
	
    fi
    
}
	       
