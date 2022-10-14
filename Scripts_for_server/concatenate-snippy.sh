#!/bin/bash


if [ $# -eq 0 ]; then
	echo "Usage concatenate-snippy.sh [options] <snippy-dir1> <snippy-dir2> ..."
	echo "Usage:"
	echo "  -f Remove duplicate mutations (default: false)"
	exit
fi

# parse option -f
while getopts 'f' option; do
        case $option in
                f)
                        filter="true"
                        ;;
        esac
done

# skip over the processed option to pass arguments
shift $((OPTIND-1))

#.......................................................

# remove existing all_snps* files

rm all_snps*

# get heading from the first file

head -n 1 $1/*tab | sed 's/^/ISOLATE\t/' >> all_snps.tab

for dir in $@; do
#    echo "Concatenate snippy folder $dir..."
	iso=$(basename $dir)
	file=$(readlink -e $dir/*tab)
	tail -n +2 $file | sed "s/^/$iso\t/"
#	(head -n 1 $file | sed 's/^/ISOLATE\t/' ; tail -n +2 $file | sed "s/^/$iso\t/")
done >> all_snps.tab



# optional filtering
if [[ $filter == "true" ]]; then
    Rscript ~/R_functions/remove_mutations_script.R all_snps.tab 


fi

 


