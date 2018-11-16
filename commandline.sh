#!/bin/bash

#---------------------------------------#
#              FUNCTIONS                #
#---------------------------------------#

function error_exit
{
  echo -e "$1" 1>&1
  exit 1
}

#---------------------------------------#
#              FILE INPUT               #
#---------------------------------------#

if [ $# -lt 14 ]
then
    error_exit "Not enough arguments provided!\nArgument list:\n[FILENAME]\n[FILE_DIR]\n[VCF_ORDER (TN or NT)]\n[Diploid_Chr]\n[Median]\n[Min_Scale]\n[Max_Scale]\n[Minx_Scale]\n[Maxx_Scale]\n[Working_Directory]\n[R_Loc]\n[RScript_Loc]\n[StartChr]\n[EndChr]\n aborting!"
fi

FILENAME=$1
FILE_DIR=$2
VCF_ORDER=$3
DIPLOID_CHR=$4
MEDIAN=$5
MINSF=$6
MAXSF=$7
XMINSF=$8
XMAXSF=$9
WORK_DIR=${10}
R_LOC=${11}
RSCRIPT_LOC=${12}
START_CHR=${13}
END_CHR=${14}

echo "Set Base Directory"

BASE_DIR=$PWD
SAMPNAME=${FILENAME%%.*}

if [ ! -d "$WORK_DIR" ]; then
	mkdir -p $WORK_DIR
fi

echo "Initialize input files"

GOOD_BAD="$BASE_DIR/vcf2cna_prep/good.bad.new"
WINDOW="$BASE_DIR/vcf2cna_prep/hg19_winbin_100bp.txt"
CHR_SIZES="$BASE_DIR/vcf2cna_prep/chr_sizes_hg19.txt"
ABB_CHR_SIZES="$BASE_DIR/vcf2cna_prep/chr.sizes.txt"
CONSPREP="$BASE_DIR/vcf2cna_prep/consprep"
SNVCOUNTS="$BASE_DIR/vcf2cna_prep/snvcounts"
HEADER="$WORK_DIR/header.txt"
BEDGRAPH="$BASE_DIR/vcf2cna_prep/bedGraphToBigWig"

# determine input filetype
echo "Determine input filetype"
head -n 1 $FILE_DIR/$FILENAME > $HEADER
FILETYPE=$(perl ${BASE_DIR}/source/file_type.pl $HEADER)
echo "$FILETYPE"

if [[ "$FILETYPE" == "unknown" ]];
then
    error_exit "Unrecognized filetype. Aborting!"
fi

# parse files based on filetype and generate snvcounts_outputfile with median_outputfile
echo "Parse files based on filetype and generate snvcounts_outputfile with median_outputfile"
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
    echo "Starting snvcounts"
    $SNVCOUNTS $FILE_DIR/$FILENAME $WORK_DIR/snvcounts_outputfile $WORK_DIR/median_outputfile
fi

# consprep
echo "Starting consprep"

if [ "$MEDIAN" == "-1" ]
then
    MEDIAN=`cat $WORK_DIR/median_outputfile`
fi

# Run CONSPREP Program and catch errors
if $CONSPREP -median=$MEDIAN -minfactor=$MINSF -maxfactor=$MAXSF -xminfactor=$XMINSF -xmaxfactor=$XMAXSF $GOOD_BAD $WINDOW $WORK_DIR/$FILENAME < $WORK_DIR/snvcounts_outputfile; then
    echo "Successfully ran consprep"
else
    error_exit "consprep crashed! aborting."
fi


# vcf2cna
echo "running VCF2CNA script"

if [ "$DIPLOID_CHR" == "-1" ] && [ "$START_CHR" == "-1" ]  && [ "$END_CHR" == "-1" ]; then
    arguments="base_dir=\"$BASE_DIR\" SAMPLE=\"$FILENAME\" forced=T hg18=F working_directory=\"$WORK_DIR\"";
fi

if [ "$DIPLOID_CHR" -ne -1 ] && [ "$START_CHR" == "-1" ] && [ "$END_CHR" == "-1" ]; then
    arguments="base_dir=\"$BASE_DIR\" SAMPLE=\"$FILENAME\" forced=T hg18=F norm.chr=$DIPLOID_CHR working_directory=\"$WORK_DIR\"";
fi

if [ "$DIPLOID_CHR" -ne -1 ] && [ "$START_CHR" -ne -1 ] && [ "$END_CHR" -ne -1 ]; then
    arguments="base_dir=\"$BASE_DIR\" SAMPLE=\"$FILENAME\" forced=T hg18=F norm.seg=c($DIPLOID_CHR,$START_CHR,$END_CHR) working_directory=\"$WORK_DIR\"";
fi

echo "$arguments"

if ${R_LOC}R --vanilla --slave --args $arguments < ${BASE_DIR}/source/VCF2CNA.R; then
	echo "Successfully ran segment analysis"
else
	error_exit "VCF2CNA.R crashed! aborting."
fi

echo "computing b-allele frequencies"
if perl ${BASE_DIR}/source/snv2baf.pl ${WORK_DIR}/snvcounts_outputfile > $WORK_DIR/$SAMPNAME.germ.baf; then
	echo "Successfully computed b-allele frequencies"
else
	error_exit "BAF.pl crashed! aborting."
fi

echo "computing tumor purity"
arguments="base_dir=\"$BASE_DIR\" working_directory=\"$WORK_DIR\" sample.name=\"$SAMPNAME\" file.name=\"$FILENAME\" analysis.type=\"fast\"";

if ${R_LOC}R --vanilla -- slave --args $arguments < ${BASE_DIR}/source/Clarity.v2.R > ${WORK_DIR}/purity.log; then
	echo "Successfully ran purity analysis"
else
	error_exit " Clarity.v2.R crashed! aborting."
fi

PURITY=`cat $WORK_DIR/Result/purity_$SAMPNAME.txt`

echo "Tumor Purity: $PURITY"

echo "running ai_plot_cnv.r"

SNVCOUNTS_OUT="$WORK_DIR/snvcounts_outputfile"

MAP1="$WORK_DIR/Result/$FILENAME"
MAP2="_CONSERTING_Mapability_100.txt"
MAPABILITY=$MAP1$MAP2

LOH1="$WORK_DIR/Result/$FILENAME"
LOH2="_LOH_RegTree.txt"
LOH=$LOH1$LOH2

IMAGE_OUTPUT="$WORK_DIR/Result/$SAMPNAME.jpg"

if ${RSCRIPT_LOC}Rscript ${BASE_DIR}/source/ai_plot_cnv.r $CHR_SIZES $SNVCOUNTS_OUT $MAPABILITY $LOH $IMAGE_OUTPUT; then
	echo "Successfully ran plotting algorithm"
else
	error_exit "ai_plot_cnv.r crashed! aborting."
fi
echo "Process Complete"
exit
