#################################################################################
 INTRO
#################################################################################
#
# This README describes how to generate a BED file:
#  * for the Agilent SureSelect Inherited Disease V1 capturing kit.
#  * which can be used by the GCC NGS pipeline.
#
# Overview:
# 
# 1.  Software and environment
# 2.  Create .genome file for genome build
# 3.  Create BED file

#################################################################################
# 1. Software and environment
#################################################################################
#
# Set environment variables for commands:
#
export ROOT_BED_FILE_DIR="/gcc/sources/BED_Files/"
export AGILENT_SS_KIT_DIR="${ROOT_BED_FILE_DIR}/Agilent_SureSelect_Inherited_Disease_V1_S0684402/"
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
# 3.  Create BED file
#################################################################################

#
# Inherited Disease V1:
#
# Downloaded S0684402.zip from Agilent's SureDesign site @
#   https://earray.chem.agilent.com/suredesign/ 
# containing the BED files for the kit (ELID S0684402):
#   S0684402_Covered.bed      (Regions covered by probes)
#   S0684402_Padded.bed       (*_Coverage.bed flanked with 100 bp)
#
# TODO: Get fragment_targets file from Agilent
#
# NOTE1: *_fragment_targets.bed has only partial overlap with *_Covered.bed:
#        probes may be shifted slightly compared to targets to compensate for repeats and/or too high/low GC%, etc.
# NOTE2: *_fragment_targets.bed is a subset of *_Padded.bed
#
mkdir -p ${AGILENT_SS_KIT_DIR}/originals/
cd ${AGILENT_SS_KIT_DIR}/originals/
# Copy S0684402.zip to ${AGILENT_SS_KIT_DIR}/originals/
unzip S0684402.zip
mv S0684402/* ./; rmdir S0684402

#
# BED file was manually processed/filtered/checked and 20 flanking bases were added by diagnotics dept.
# Received HPO_+en-20_target_v2.bed from Martine Veldhuis <m.t.veldhuis@umcg.nl>
#
mkdir -p ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/
cd ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/
cp ${ROOT_BED_FILE_DIR}/human_g1k_v37.genome ./
export FLANKED_BED_FILE_PREFIX=S0684402_UMCG
# Copy HPO_+en-20_target_v2.bed to ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/ dir.

#
# Merge overlapping regions.
#
sort -k1,1 -k2,2n \
   ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/HPO_+en-20_target_v2.bed \
 | mergeBed \
   -nms \
   -i stdin \
 > ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed

#
# Flatten annotation.
#
perl -F'\t' -lane \
   '$F[3] = join(",", sort(grep(!$_{$_}++, split(",", $F[3])))); 
    print join "\t", @F; %_ = ();
   ' \
   ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed \
 > ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed.tmp
mv ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed{.tmp,}

#
# Count total amount of bases:
#
for bed_file in {HPO_+en-20_target_v2.bed,${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed}; do \
   awk '{total_bases+=$3-$2} END {printf "   %s", total_bases}' ${bed_file}; \
   echo " ${bed_file}";
done
#
#   8348547 HPO_+en-20_target_v2.bed
#   8343358 S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.bed
#
# Redundancy removed = 5189 bases.
#

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
mv ${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed{,.tmp}
for char_class in {'0-9','X','Y','M'}; do \
  awk "\$1 ~ /^[${char_class}]/" \
     ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed.tmp \
  |  sort -k1n -k2n \
  >> ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed; \
done
rm ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed.tmp

#
# Rename final file for consistency and readability.
#
mv ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/${FLANKED_BED_FILE_PREFIX}.sorted.stripped.flanked20bp.regionMerged.bed \
   ${AGILENT_SS_KIT_DIR}/gcc_version_indexed_for_human_g1k_v37/Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.labelled.bed

###########################################################################################################
# Final BED file used in experiments:
# Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.labelled.bed
###########################################################################################################





