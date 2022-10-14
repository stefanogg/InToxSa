# Assess arguments
if [ $# -lt  3 ]
then
	echo "Insufficient number of arguments: please provide <PAIR ID> <iso1> <iso2>"
	exit 1
fi

# Get variables
PAIR_ID=$1
iso1=$2
iso2=$3

# Starting statements
echo "Running $0 on $PAIR_ID"
echo "Reference (iso1) is $iso1"
echo "iso2 is $iso2"

# Create pair directory
if [ -d $PAIR_ID ]
then
	echo "$PAIR_ID exists. $0 will remove $PAIR_ID"
	echo "In the future we will upgrade $0 to include a --force option"
	rm -r $PAIR_ID
fi

mkdir $PAIR_ID
cd $PAIR_ID

 # Create isolates.tab file
~/perl5/bin/readprint3.sh -d ~/dropbox/Sa_ANZCOSS/ $iso1 $iso2 > isolates.tab

# Copy reference
mkdir reference 
cp ~/VANANZ_phenotypes/assemblies/prokka/$iso1/*.gbk reference

# Run snippy
cat isolates.tab | while read f1 f2 f3;
	do echo "Running snippy on $f1"
	snippy --ref reference/*.gbk --outdir $f1 --R1 $f2 --R2 $f3
done

cd ..

# Final statements
echo "$0 on $PAIR_ID completed"
echo "You can run-concatenate-snippy.sh to remove variants from the reference reads"
echo "Or, even better, you can run mask-variants.sh to remove variants from the reference reads and variants in positions where reference reads have low coverage"
echo "This will substantially improve the accuracy of your mutation list!"
