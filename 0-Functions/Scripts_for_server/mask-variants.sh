#!/usr/bin/env bash

# add options using getopt (later)

# Author: Stefano Giulieri
# Date: February 2021 (revision)

USAGE="$0 <grouped snippy dir>\\n
\\t Grouped snippy dir: a directory with 2 or more snippy runs of related sequences and a shared reference"


if [ $# -eq 0 ]
then
	echo -e $USAGE
	exit 1
fi

echo "This is mask-variants.sh" 
echo "This scripts uses bedtools and Torsten's snippy-vcf_tp_tab script to mask variants from reference reads and variants in regions where reference reads have low coverage" 

d=$(readlink -e $1) # absolute path of the grouped snippy directory
p=$(basename $d) # id of the snippy group
l="mask-variants.log"

# check if group directory exists and remove if true (later to be replaced by --force option)
if [ -d $p ]
then
	echo "Directory $p exists. Removing"
	rm -r $p
fi

mkdir -p $p # create group directory
cd $p

echo "Running mask-variants.sh on snippy group $p" 2>&1 | tee -a $l

# copy relevant reference files
ref=$(basename -s .gbk $d/reference/*.gbk)
echo "Copying reference $ref files snps.bed snps.aligned.fa.bed ref.fa ref.gff" 2>&1 | tee -a $l

 
cp $d/$ref/snps.bed ref.snps.bed
cp $d/$ref/snps.aligned.fa.bed ref.snps.aligned.fa.bed
#cp $d/$ref/ref.fa . # shortcut if these files are present
#cp $d/$ref/ref.gff .
cp $d/reference/*.gbk ref.gbk
any2fasta -u ref.gbk > ref.fa
genbank2gff.pl ref.gbk > ref.gff

# mask variants
cp $d/isolates.tab .

for i in $(cut -f1 isolates.tab)
do
	echo "Masking variants for strain $i"
	snippy_d=$(readlink -e $d/$i)
	mkdir -p $i
	cd $i
	cp $snippy_d/snps.vcf .
	cp $snippy_d/snps.tab .
	bedtools subtract -A -a snps.vcf -b ../ref.snps.bed -header > snps.noref.vcf # we need to use the option -A to remove longer mutations (eg complex, indels)
	bedtools subtract -a snps.noref.vcf -b ../ref.snps.aligned.fa.bed -header > snps.mask.vcf

	snippy-vcf_to_tab --ref ../ref.fa --gff ../ref.gff -vcf snps.noref.vcf > snps.noref.tab	
	snippy-vcf_to_tab --ref ../ref.fa --gff ../ref.gff -vcf snps.mask.vcf > snps.mask.tab

	cd ..
done 2>&1 | tee -a $l

# concatenate variants in the group

HEADER="CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE\tFTYPE\tSTRAND\tNT_POS\tAA_POS\tEFFECT\tLOCUS_TAG\tGENE\tPRODUCT"
REFGBK=$ref.gbk

for i in snps.tab snps.noref.tab snps.mask.tab
do
	echo -e "PAIR_ID\tREFERENCE\tISOLATE\t$HEADER" >> all_$i
	for j in $(cut -f1 isolates.tab)
	do
		sed '1d' $j/$i | sed "s/^/$p\t$REFGBK\t$j\t/" 
	done >> all_$i
done

# trick to add columns using csvtk mutate2
# val=ABC && seq 3 | csvtk mutate2 -H -e " '$val' "

