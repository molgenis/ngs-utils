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

#################################################################################
# 1. Software and environment
#################################################################################
#
# Set environment variables for commands:
#
export ROOT_BED_FILE_DIR="/gcc/sources/BED_Files/"
export ASS5_DIR="${ROOT_BED_FILE_DIR}/Agilent_SureSelect_Human_All_Exon_V5_S04380110/"
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
# Downloaded S04380110.zip from Agilent's SureDesign site @
#   https://earray.chem.agilent.com/suredesign/ 
# containing the BED files for the All Exon version 5 kit (ELID S04380110):
#   S04380110_Covered.bed      (Regions covered by probes)
#   S04380110_Padded.bed       (*_Coverage.bed flanked with 100 bp)
#
# Received via email from Pierre Bourbon, PhD – Sequencing Applications Support Specialist <pbourbon@agilent.com>:
#   SS_V5_fragment_targets.bed (Target regions selected to design the probes)
#
# NOTE1: *_fragment_targets.bed has only partial overlap with *_Covered.bed:
#        probes may be shifted slightly compared to targets to compensate for repeats and/or too high/low GC%, etc.
# NOTE2: *_fragment_targets.bed is a subset of *_Padded.bed
#
export EXOME_PADDED_BED_FILE_PREFIX=S04380110_Padded
export EXOME_TARGET_BED_FILE_PREFIX=SS_V5_fragment_targets

mkdir -p ${ASS5_DIR}/originals/
cd ${ASS5_DIR}/originals/
# Copy S04380110.zip to ${ROOT_SS5_DIR}/originals/
unzip S04380110.zip
mv S04380110/* ./; rmdir S04380110
# Copy SS_V5_fragment_targets.bed to ${ROOT_SS5_DIR}/originals/ dir.
mkdir -p ${ASS5_DIR}/gcc_version_indexed_for_human_g1k_v37/
cd ${ASS5_DIR}/gcc_version_indexed_for_human_g1k_v37/
cp ${ROOT_BED_FILE_DIR}/human_g1k_v37.genome ./

#
# Remove chr prefix for chromosome names, remove leading "track" + "browser" lines and sort
#
sed 's/^chr//' ${ASS5_DIR}/originals/${EXOME_PADDED_BED_FILE_PREFIX}.bed \
 | sed 1,2d \
 | sort -k1,1 -k2,2n \
 > ${ASS5_DIR}/gcc_version_indexed_for_human_g1k_v37/${EXOME_PADDED_BED_FILE_PREFIX}.sorted.bed

sed 's/^chr//' ${ASS5_DIR}/originals/${EXOME_TARGET_BED_FILE_PREFIX}.bed \
 | sort -k1,1 -k2,2n \
 > ${ASS5_DIR}/gcc_version_indexed_for_human_g1k_v37/${EXOME_TARGET_BED_FILE_PREFIX}.sorted.bed

#
# Add 20 flanking bases to Targets
#
slopBed -b 20 \
-i ${EXOME_TARGET_BED_FILE_PREFIX}.sorted.bed \
-g human_g1k_v37.genome \
| mergeBed \
  -i stdin \
 > ${EXOME_TARGET_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed

#
# Drop annotation from PADDED, because 
#   * it is messy, 
#   * the TARGET BED file does not have any annotation and
#   * mergeBed will fail to merge BED files when the amount of columns is not the same.
# and merge (partially) overlapping regions using mergeBed
#   * NOTE: always sort sequence name and position with sort -k1,1 -k2,2n before merging to prevent mistakes in merging!
#
awk -F "\t" 'BEGIN {OFS=FS} {print $1,$2,$3}' ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.bed \
| sort -k1,1 -k2,2n \
| mergeBed \
  -i stdin \
  > ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.regionMerged.bed

#
# Merge TARGET with PADDED.
#   * NOTE: always sort sequence name and position with sort -k1,1 -k2,2n before merging to prevent mistakes in merging!
#
cat ${EXOME_TARGET_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed \
    ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.regionMerged.bed \
| sort -k1,1 -k2,2n \
| mergeBed \
  -i stdin \
  >  ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.bed

#
# Count total amount of bases to make sure we added coverage and did not loose any:
#
for bed_file in {${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.regionMerged.bed,${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.bed}; do \
   awk '{total_bases+=$3-$2} END {printf "   %s", total_bases}' ${bed_file}; \
   echo " ${bed_file}";
done
#
#   89476292 S04380110_Padded.sorted.stripped.regionMerged.bed
#   89567011 S04380110_Padded.sorted.stripped.flanked0-20bp.regionMerged.bed
#
#   Added bases = (89567011-89476292)/89476292*100 = 0.1%
#

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
mv ${EXOME_PADDED_BED_FILE_PREFIX}.{sorted.,}stripped.flanked0-20bp.regionMerged.bed
awk '$1 ~ /^[0-9]/' ${EXOME_PADDED_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.bed | sort -k1n -k2n >  ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.bed
awk '$1 ~ /^[X]/'   ${EXOME_PADDED_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.bed | sort -k1n -k2n >> ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.bed
awk '$1 ~ /^[Y]/'   ${EXOME_PADDED_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.bed | sort -k1n -k2n >> ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.bed
awk '$1 ~ /^[M]/'   ${EXOME_PADDED_BED_FILE_PREFIX}.stripped.flanked0-20bp.regionMerged.bed | sort -k1n -k2n >> ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.bed

#
# Label regions.
#
awk -F "\t" 'BEGIN {OFS=FS} {print $1,$2,$3,"Agilent_SS_AllExonV5"}' \
   ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.bed \
 > ${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed

#################################################################################
# Used in Exome Seq experiments for:	File:
#################################################################################
# 1. Variant calling + Coverage BEDs	${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
# 2. QC metrics Targets					${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
# 3. QC metrics capt. kit performance	${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.regionMerged.bed
#################################################################################





