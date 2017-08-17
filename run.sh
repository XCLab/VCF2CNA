#!/bin/bash

#---------------------------------------#
#              FUNCTIONS                #
#---------------------------------------#

function error_exit
{
  echo "$1" 1>&2
  exit 1
}

#---------------------------------------#
#              FILE INPUT               #
#---------------------------------------#

if [ $# -eq 4 ]
then
	FILENAME=$1
	FILE_DIR=$2
	SEQ_TYPE=$3
	VCF_ORDER=$4
fi

if [ $# -eq 3 ]
then
	FILENAME=$1
	FILE_DIR=$2
	SEQ_TYPE=$3
	VCF_ORDER="NA"
fi

if [ $# -lt 3 ]
then
	error_exit "[FILENAME][FILE_DIR][SEQ_TYPE (EXOME or WHOLEGENOME)]! aborting."
fi

FILENAME=$1
FILE_DIR=$2
SEQ_TYPE=$3

if [[ "$SEQ_TYPE" != "EXOME" && "$SEQ_TYPE" != "WHOLEGENOME" ]];
then
	error_exit "Incorrect seq_type: ${SEQ_TYPE}.  Options are EXOME or WHOLEGENOME."
fi

BASE_DIR=$PWD
SAMPNAME=${FILENAME%%.*}
WORK_DIR="$BASE_DIR/$SAMPNAME/$SEQ_TYPE"

if [ ! -d "$WORK_DIR" ]; then
	mkdir -p $WORK_DIR
fi

module load R

GOOD_BAD="$BASE_DIR/vcf2cna_prep/good.bad.new"
WINDOW="$BASE_DIR/vcf2cna_prep/hg19_winbin_100bp.txt"
CHR_SIZES="$BASE_DIR/vcf2cna_prep/chr_sizes_hg19.txt"
CONSPREP="$BASE_DIR/vcf2cna_prep/consprep"
SNVCOUNTS="$BASE_DIR/vcf2cna_prep/snvcounts"
HEADER="$WORK_DIR/header.txt"

# determine input filetype
head -n 1 $FILE_DIR/$FILENAME > $HEADER
FILETYPE=$(perl ${BASE_DIR}/source/file_type.pl $HEADER)
echo "$FILETYPE"

# parse files based on filetype and generate snvcounts_outputfile with median_outputfile
if [ "$FILETYPE" == "VCF" ]
then
	if [ $# -lt 3 ]
	then
		error_exit "For VCF_FILES you must specify pair order: TN or NT! aborting."
	fi
	VCF_SNPCOUNT=$(perl ${BASE_DIR}/source/vcf_parser_4.1.pl $FILE_DIR/$FILENAME $VCF_ORDER $WORK_DIR)
	grep 'Pos' $WORK_DIR/snvcounts_outputfile > $WORK_DIR/sort_file
	grep -v -E 'X|Y|Pos' $WORK_DIR/snvcounts_outputfile | sort -n -k1.4 -k2 >> $WORK_DIR/sort_file
	grep -E 'X|Y' $WORK_DIR/snvcounts_outputfile | sort -t$'\t' -k1.4,1d -k2,2n >> $WORK_DIR/sort_file
	mv $WORK_DIR/sort_file $WORK_DIR/snvcounts_outputfile

fi

if [ "$FILETYPE" == "SJ_MAF" ]
then
	cp $FILE_DIR/$FILENAME $WORK_DIR/snvcounts_outputfile
	HEADER_LINE="Chr\tPos\tTumorMutant\tTumorTotal\tNormalMutant\tNormalTotal"
	sed -i "1s/.*/$HEADER_LINE/" $WORK_DIR/snvcounts_outputfile
	SJ_MAF_MEDIAN=$(perl ${BASE_DIR}/source/sj_maf_parser.pl $WORK_DIR/snvcounts_outputfile > $WORK_DIR/median_outputfile)
fi

if [[ "$FILETYPE" == "HIGH20" || "$FILETYPE" == "MAF" ]];
then
	$SNVCOUNTS $FILE_DIR/$FILENAME $WORK_DIR/snvcounts_outputfile $WORK_DIR/median_outputfile
fi

# consprep
echo "Starting consprep"
MEDIAN=`cat $WORK_DIR/median_outputfile`
if $CONSPREP -median=$MEDIAN $GOOD_BAD $WINDOW $WORK_DIR/$FILENAME < $WORK_DIR/snvcounts_outputfile; then
	echo "Successfully ran consprep"
else
	error_exit "consprep crashed! aborting."
fi

# vcf2cna
echo "running VCF2CNA script"
arguments="SAMPLE=\"$FILENAME\" forced=T hg18=F working_directory=\"$WORK_DIR\"";

if R --vanilla --slave --args $arguments < ${BASE_DIR}/source/VCF2CNA.R; then
	echo "Successfully ran segment analysis"
else
	error_exit "VCF2CNA.R crashed! aborting."
fi

echo "running ai_plot_cnv.r"

SNVCOUNTS_OUT="$WORK_DIR/snvcounts_outputfile"

MAP1="$WORK_DIR/Result/$FILENAME"
MAP2="_CONSERTING_Mapability_100.txt"
MAPABILITY=$MAP1$MAP2

LOH1="$WORK_DIR/Result/$FILENAME"
LOH2="_LOH_RegTree.txt"
LOH=$LOH1$LOH2

IMAGE_OUTPUT="$WORK_DIR/Result/$SAMPNAME.jpg"

if Rscript ${BASE_DIR}/source/ai_plot_cnv.r $CHR_SIZES $SNVCOUNTS_OUT $MAPABILITY $LOH $IMAGE_OUTPUT; then
	echo "Successfully ran plotting algorithm"
else
	error_exit "ai_plot_cnv.r crashed! aborting."
fi

echo "computing b-allele frequencies"
if perl ${BASE_DIR}/source/snv2baf.pl ${WORK_DIR}/snvcounts_outputfile > $WORK_DIR/$SAMPNAME.germ.baf; then
	echo "Successfully computed b-allele frequencies"
else
	error_exit "BAF.pl crashed! aborting."
fi

echo "computing tumor purity"
arguments="working_directory=\"$WORK_DIR\" sample.name=\"$SAMPNAME\" file.name=\"$FILENAME\" analysis.type=\"fast\"";

if R --vanilla -- slave --args $arguments < ${BASE_DIR}/source/Clarity.v2.R > ${WORK_DIR}/purity.log; then
	echo "Successfully ran purity analysis"
else
	error_exit " Clarity.v2.R crashed! aborting."
fi


