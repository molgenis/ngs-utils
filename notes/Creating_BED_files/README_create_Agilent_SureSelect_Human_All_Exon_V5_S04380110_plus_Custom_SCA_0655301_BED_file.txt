#################################################################################
 INTRO
#################################################################################
#
# This README describes how to generate a BED file:
#  * for the Agilent SureSelect Human All Exon V5 capturing kit.
#  * which can be used by the GCC NGS pipeline.
#
# Overview:
# 
# 1.  Software and environment
# 2.  Create .genome file for genome build
# 3.  Create Agilent SureSelect All Exon V5 capturing kit BED file
# 4.  Supplement BED file with Custom SCA targets

#################################################################################
# 1. Software and environment
#################################################################################
#
# Set environment variables for commands:
#
export ROOT_BED_FILE_DIR="/gcc/sources/BED_Files/"
export ASS5_DIR="${ROOT_BED_FILE_DIR}/Agilent_SureSelect_Human_All_Exon_V5_S04380110/"
export ASS5_PLUS_SCA_DIR="${ROOT_BED_FILE_DIR}/Agilent_SureSelect_Human_All_Exon_V5_S04380110_plus_Custom_SCA_0655301/"
# 
# The following software was installed on cluster.gcc.rug.nl
#  * bedtools2-2.19.1 from https://github.com/arq5x/bedtools2/releases
# and added to the environment using:
#
module load bedtools/2.19.1

#################################################################################
# 2. Create .genome file for genome build
#################################################################################
#
# Create *.genome file from the dictionary file of our currently used reference genome
# This *.genome file contains the start and stop of all chromosomes.
# When manipulating a *.BED file with bedtools and a region falls outside the chromosome start/end, 
# it will automatically be clipped/trimmed to the chromosome boundary.

#
# Get copy of the *.dict file for our reference genome.
#
cd ${ROOT_BED_FILE_DIR}
cp /gcc/resources/b37/indices/human_g1k_v37.dict ${ROOT_BED_FILE_DIR}

#
# Skip headerline and print second + third column from reference dictionary file:
#
sed '1,1d' ${ROOT_BED_FILE_DIR}/human_g1k_v37.dict | \
awk '{print $2,$3}' OFS="\t" > \
${ROOT_BED_FILE_DIR}/human_g1k_v37.genome

#
# Replace SN and LN tags with "nothing"
#
perl -pi -e 's/SN://g' ${ROOT_BED_FILE_DIR}/human_g1k_v37.genome
perl -pi -e 's/LN://g' ${ROOT_BED_FILE_DIR}/human_g1k_v37.genome

#################################################################################
# 3.  Create Agilent SureSelect All Exon V5 capturing kit BED file
#################################################################################

#
# Exome V5:
#
# Construction of BED file for the All Exon version 5 kit (ELID S04380110) is described in detail in:
# README_create_Agilent_SureSelect_Human_All_Exon_V5_S04380110_BED_file.txt
#
export EXOME_PADDED_BED_FILE_PREFIX=S04380110_Padded

#################################################################################
# 4.  Supplement BED file with Custom SCA targets
#################################################################################

#
# Custom probes for SCA targets were designed with:
#   https://earray.chem.agilent.com/suredesign/ 
# resulting in Design/kit ELID 0655301.
#
# Received BED file for 0655301 via email from Esther Nibbeling <e.a.r.nibbeling@umcg.nl>:
#   SCA3 plus targets.bed
# Targets do not need additional flanking bases.
# 
export SCA_TARGET_BED_FILE_PREFIX="Agilent_SureSelect_Human_Custom_SCA_0655301"
export EXOME_PLUS_SCA_BED_FILE_PREFIX="Agilent_SureSelect_Human_All_Exon_V5_S04380110_plus_Custom_SCA_0655301"

mkdir -p ${ASS5_PLUS_SCA_DIR}/originals/
cd ${ASS5_PLUS_SCA_DIR}/originals/
# Copy SCA3 plus targets.bed to ${ASS5_PLUS_SCA_DIR}/originals/
cp 'SCA3 plus targets.bed' ${ASS5_PLUS_SCA_DIR}/originals/SCA3_plus_targets.bed
mkdir -p ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/
cd ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/
cp ${ROOT_BED_FILE_DIR}/human_g1k_v37.genome ./

#
# Remove chr prefix for chromosome names, drop annotation column 4 and sort
#
sed 's/^chr//' ${ASS5_PLUS_SCA_DIR}/originals/SCA3_plus_targets.bed \
 | cut -f 1,2,3 \
 | sort -k1,1 -k2,2n \
 > ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/${SCA_TARGET_BED_FILE_PREFIX}.sorted.bed

#
# Label regions.
#
awk -F "\t" 'BEGIN {OFS=FS} {print $1,$2,$3,"Agilent_SS_SCA3"}' \
   ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/${SCA_TARGET_BED_FILE_PREFIX}.sorted.bed \
 > ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/${SCA_TARGET_BED_FILE_PREFIX}.sorted.labelled.bed


#
# Merge standard Exome V5 with custom SCA targets.
#   * NOTE: always sort sequence name and position with sort -k1,1 -k2,2n before merging to prevent mistakes in merging!
#
cat ${ASS5_DIR}/gcc_version_indexed_for_human_g1k_v37/${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed \
    ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/${SCA_TARGET_BED_FILE_PREFIX}.sorted.labelled.bed \
| sort -k1,1 -k2,2n \
| mergeBed \
  -nms \
  -i stdin \
  >  ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed

#
# Count total amount of bases to make sure we added coverage and did not loose any:
#
declare -a BED_FILES=(
  ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/${SCA_TARGET_BED_FILE_PREFIX}.sorted.labelled.bed
  ${ASS5_DIR}/gcc_version_indexed_for_human_g1k_v37/${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
  ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
)
for bed_file in ${BED_FILES[@]}; do \
   awk '{total_bases+=$3-$2} END {printf "   %s", total_bases}' ${bed_file}; \
   echo " ${bed_file}";
done
#
#   1165889  Agilent_SureSelect_Human_Custom_SCA_0655301.sorted.labelled.bed
#   89567011 S04380110_Padded.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
#   90578949 Agilent_SureSelect_Human_All_Exon_V5_S04380110_plus_Custom_SCA_0655301.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
#
#   Added bases = (90578949-89567011)/89567011*100 = 1.13%
#

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
cd ${ASS5_PLUS_SCA_DIR}/gcc_version_indexed_for_human_g1k_v37/
mv ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.{sorted.,}stripped.flanked0-20bp.regionMerged.labelled.bed
awk '$1 ~ /^[0-9]/' ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.labelled.bed | sort -k1n -k2n >  ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
awk '$1 ~ /^[X]/'   ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.labelled.bed | sort -k1n -k2n >> ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
awk '$1 ~ /^[Y]/'   ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.labelled.bed | sort -k1n -k2n >> ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
awk '$1 ~ /^[M]/'   ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.labelled.bed | sort -k1n -k2n >> ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
rm ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.labelled.bed

#
# Flatten annotation.
#
# Add 'ASD_Symbol=' prefix and
# remove duplicate values in one go using perl:
#
perl -F'\t' -lane \
   '$F[3] = join(",", sort(grep(!$_{$_}++, split(",", $F[3])))); 
    print join "\t", @F; %_ = ();
   ' \
   ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed \
 > ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed.flattened \
 ; mv ${EXOME_PLUS_SCA_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed{.flattened,}
 
###########################################################################################################
# Final BED file used in Exome Seq experiments:
# Agilent_SureSelect_Human_All_Exon_V5_S04380110_plus_Custom_SCA_0655301.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
###########################################################################################################





