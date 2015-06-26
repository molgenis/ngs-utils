#################################################################################
 INTRO
#################################################################################
#
# This README describes how to generate a BED file:
#  * for genes based on official HGNC Gene Symbols; 
#  * in this case for genes from the Clinical Genomics Database (CGD)
#  * which can be used by the GCC NGS pipeline.
#
# Overview:
# 
# 1.  Software and environment
# 2.  Create .genome file for genome build
# 3.A Create CGD BED file using custom Perl script and Ensembl API (Ensembl)
# 3.B Create CGD BED file using Agilent SureDesign (ASD) tool
# 4.  Compare Ensembl and ASD BED files and merge into single CGD BED file
# 5.  Compare CGD and Agilent SureSelect All Exon V5 capturing kit BED files 
#     and merge into single Agilent SureSelect All Exon V5 plus CGD BED file.

#################################################################################
# 1.  Software and environment
#################################################################################
#
# Set environment variables for commands:
#
export ROOT_BED_FILE_DIR="/gcc/sources/BED_Files/"
export CGD_DATE="2014-10-08"
export CGD_FILE="CGD_${CGD_DATE}_updated.txt"
export CGD_DIR="${ROOT_BED_FILE_DIR}/ClinicalGenomicsDatabase_${CGD_DATE}/"
export ENSEMBL_VERSION=75
#
# Choose transcript types to target:
#  * all: regions of all transcripts of the selected genes are included in the BED file.
#  * CT:  complete regions of protein Coding Transcripts are included. Hence UTRs are included.
#  * CDS: only the the Coding DNA Sequences of protein coding transcripts are included. Hence UTRs are excluded.
# 
# We chose for all: variants in any non-protein coding parts of genes, 
# which will be hard to interpret will in most cases be filtered out as "VOUS" 
# when filtering in downstream software like Cartagenia.
#
export TARGET="all"
#export TARGET="CDS"
#export TARGET="CT"
# 
# The following software was used
#  * bedtools2-2.19.1 from https://github.com/arq5x/bedtools2/releases
#  * Ensembl API version 75
#  * bioperl-1.2.3 (required for Ensembl API)
#  * http://www.bbmriwiki.nl/svn/ngs_scripts/trunk/GeneSymbol2BED.pl
#  * Online Agilent SureDesign tool @ https://earray.chem.agilent.com/suredesign/
#
# and added to the environment using:
#
module load bedtools/2.19.1

#################################################################################
# 2.  Create .genome file for genome build
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
# 3.A Create CGD BED file using custom Perl script and Ensembl API (Ensembl)
#################################################################################

#
# Download a copy of the CGD and update/check the gene symbols.
# For details refer to:
#    README_CGD_${CGD_DATE}_update_GeneSymbols.txt
# in the same location as this README.
# This should result in a ${CGD_FILE}
#

#
# Create subdir.
#
mkdir -p ${CGD_DIR}/BED/From_Ensembl_API/
cd       ${CGD_DIR}/BED/From_Ensembl_API/

#
# Get Perl script to query ENSEMBL
#
#  * ENSEMBL Perl API needs to be installed first (this depends on Bioperl)
#  * Input must be TAB delimited file with header line and at least these 2 columns: 
#      * HGNC Gene Symbols  -  default column label: "GENE"
#      * HGNC Gene IDs      -  default column label: "HGNC ID"
#    Columns will be detected based on column label in a header line.
#    Additional columns may be present.
#    When columns have a different label these can be specified with commandline arguments.
#
curl http://www.bbmriwiki.nl/svn/ngs_scripts/trunk/GeneSymbol2BED.pl > GeneSymbol2BED.pl

#
# Optionally for debugging you can try a small subset:
#
#export CGD_FILE_SUBSET="CGD_10_genes_2014-10-08_updated.txt"
#head -11 ${CGD_DIR}/${CGD_FILE} > ${CGD_DIR}/${CGD_FILE_SUBSET}
#
#./GeneSymbol2BED.pl \
#-i ${CGD_FILE_SUBSET} \
#-g HGNC_Symbol \
#-h HGNC_ID \
#-o ${CGD_FILE_SUBSET}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.bed \
#-t ${TARGET}

#
# Final run for the complete CGD:
#
./GeneSymbol2BED.pl \
-i ${CGD_DIR}/${CGD_FILE} \
-g 'Approved Symbol' \
-h 'HGNC ID' \
-o ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.bed \
-t ${TARGET} \
 > ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.log 2>&1
#
# Log:
#
# 2014/10/27 23:50:58 L:184 INFO> Parsing HGNC Gene Symbols from CGD_2014-10-08_updated.txt...
# 2014/10/27 23:50:58 L:185 INFO> and storing (genomic) regions for targets "all" to CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75.bed...
# 2014/10/28 00:00:14 L:332 ERROR> Cannot find gene FSHMD1A in Ensembl core database!
# 2014/10/28 00:00:21 L:332 ERROR> Cannot find gene ACKR1 in Ensembl core database!
# 2014/10/28 00:16:48 L:332 ERROR> Cannot find gene BCO1 in Ensembl core database!
# 2014/10/28 00:18:13 L:332 ERROR> Cannot find gene KIZ in Ensembl core database!
# 2014/10/28 00:20:09 L:332 ERROR> Cannot find gene CEP83 in Ensembl core database!
# 2014/10/28 00:25:19 L:332 ERROR> Cannot find gene CFAP57 in Ensembl core database!
# 2014/10/28 00:25:22 L:332 ERROR> Cannot find gene CFAP53 in Ensembl core database!
# 2014/10/28 00:27:48 L:332 ERROR> Cannot find gene ATXN8 in Ensembl core database!
# 2014/10/28 00:27:54 L:332 ERROR> Cannot find gene KCNJ18 in Ensembl core database!
# 2014/10/28 00:27:55 L:113 INFO> ===================================================================
# 2014/10/28 00:27:55 L:114 INFO> Genes Symbols found and represented in BED file: 2924.
# 2014/10/28 00:27:55 L:115 INFO> Genes Symbols discarded and lacking in BED file: 0. (Discarded because sequence not part of target all).
# 2014/10/28 00:27:55 L:116 INFO> Genes Symbols missing and lacking in BED file:   9.
# 2014/10/28 00:27:55 L:117 INFO> ===================================================================
# 2014/10/28 00:27:55 L:118 INFO> Finished!

#
# For some of the problematic genes we can include their regions using a Previous Symbol or Alias:
# 
# HGNC     Previous Symbol or Alias present in Ensembl75
# ACKR1    DARC
# BCO1     BCMO1
# KIZ      PLK1S1
# CEP83    CCDC41
# CFAP57   WDR65
# CFAP53   CCDC11
# ATXN8    ---
# FSHMD1A  ---
# KCNJ18   ---
#
head -1 ${CGD_DIR}/${CGD_FILE} > ${CGD_DIR}/${CGD_FILE}.subset_with_alias.txt
for genesymbol in {'ACKR1','BCO1','KIZ','CEP83','CFAP53','CFAP57'}; do fgrep ${genesymbol} ${CGD_DIR}/${CGD_FILE} >> ${CGD_DIR}/${CGD_FILE}.subset_with_alias.txt; done
# Add ALIAS column manually to ${CGD_DIR}/${CGD_FILE}.subset_with_alias.txt and then
./GeneSymbol2BED.pl \
-i ${CGD_DIR}/${CGD_FILE}.subset_with_alias.txt \
-g 'ALIAS' \
-h 'HGNC ID' \
-o ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.bed.extra \
-t ${TARGET} \
 > ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.log.extra 2>&1
#
# Log:
#
2014/11/07 14:12:52 L:184 INFO> Parsing HGNC Gene Symbols from CGD_2014-10-08_updated.txt.subset_with_alias.txt...
2014/11/07 14:12:52 L:185 INFO> and storing (genomic) regions for targets "all" to CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75.bed.extra...
2014/11/07 14:12:59 L:113 INFO> ===================================================================
2014/11/07 14:12:59 L:114 INFO> Genes Symbols found and represented in BED file: 6.
2014/11/07 14:12:59 L:115 INFO> Genes Symbols discarded and lacking in BED file: 0. (Discarded because sequence not part of target all).
2014/11/07 14:12:59 L:116 INFO> Genes Symbols missing and lacking in BED file:   0.
2014/11/07 14:12:59 L:117 INFO> ===================================================================
2014/11/07 14:12:59 L:118 INFO> Finished!
#
# Append regions for genes included by a gene symbol alias.
#
cat ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.bed.extra \
 >> ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.bed

#
# Flank, sort and merge.
#
# Flank regions in BED file with 20bp on both ends using bedtools
# and use generated *.genome file to trim regions to chromosomal boundaries where necessary.
# Note: must flank first before region merge otherwise new overlapping regions may appear after flanking!
#
slopBed -b 20 \
  -i ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.bed \
  -g ${ROOT_BED_FILE_DIR}/human_g1k_v37.genome \
| sort -k1,1V -k2,2n \
| mergeBed -nms \
  -i stdin \
>  ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.bed

#
# Simplify annotation by flattening.
#
# First remove ENSE* exon IDs and ENST* transcript IDs using sed
# and then remove duplicate values using perl:
#
sed 's/ENSE[0-9]*|ENST[0-9]*|//g' \
  ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.bed \
| perl -F'\t' -lane'
  $F[3] = join ",", grep !$_{$_}++, split ",", $F[3];
  print join "\t", @F; %_ = ();
  ' \
> ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
mv ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed{,.tmp}
for char_class in {'0-9','X','Y','M'}; do \
  awk "\$1 ~ /^[${char_class}]/" \
  ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp \
  | sort -k1V -k2n \
  >> ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed; \
done
rm ${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp

#################################################################################
# 3.B Create CGD BED file using Agilent SureDesign (ASD) tool
#################################################################################

#
# Download a copy of the CGD and update/check the gene symbols.
# For details refer to:
#    README_CGD_${CGD_DATE}_update_GeneSymbols.txt
# in the same location as this README.
# This should result in a ${CGD_FILE}
#

#
# Create subdir.
#
mkdir -p ${CGD_DIR}/BED/From_ASD/
cd       ${CGD_DIR}/BED/From_ASD/

#
# Use HGNC Gene Symbols from ${CGD_FILE} to also create a custom capturing kit with the Agilent SureDesign tool @
#
#   https://earray.chem.agilent.com/suredesign/
#
#    * Specify the relevant regions to include: only CDS, including UTR or all transcripts
#    * DO NOT flank regions using ASD tool: ASD has all starts wrong (+1) as it fails to adjust starts for zero-based BED format.
#
# and download the BED file for this kit as ${CGD_FILE}.${TARGET}-targetsFromASD.bed
#

#
# When the ASD BED file was made on a Windows computer you may have to fix
#  1. convert UTF16 -> UTF8 encoding
#  2. convert the line end characters
# in exactly that order!
#
iconv -f UTF-16 -t UTF-8 ${CGD_FILE}.${TARGET}-targetsFromASD.bed > ${CGD_FILE}.${TARGET}-targetsFromASD.bed.utf8
mv ${CGD_FILE}.${TARGET}-targetsFromASD.bed{.utf8,}
dos2unix ${CGD_FILE}.${TARGET}-targetsFromASD.bed

#
# Remove chr prefix for chromosomes.
#
perl -pi -e 's/^chr//' ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.bed

#
# Drop patches, flank, sort and merge.
#
# Drop any sequences/chromosomes that do not start with 0-9, X, Y, or M.
# This removes regions located on extra chromosomal super contigs; usually
# these are alternative assemblies or patches.
# These regions cannot be processed with the current genome assembly as used in GCC pipelines 
# and need to be removed to prevent errors / crashing jobs.
#
# Flank regions in BED file with -20bp and +20bp using bedtools
# and use generated *.genome file to trim regions to chromosomal boundaries where necessary.
# Note: must flank first before region merge otherwise new overlapping regions may appear after flanking!
# Note: ASD used to have all starts wrong (+1) as it failed to adjust starts for zero-based BED format,
#       but this appears te be solved: no need to compensate anymore by flanking with -21bp.
#
grep -i '^[0-9XYM]' \
  ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.bed \
| slopBed -l 20 -r 20 \
  -i stdin \
  -g ${ROOT_BED_FILE_DIR}/human_g1k_v37.genome \
| sort -k1,1V -k2,2n \
| mergeBed -nms \
  -i stdin \
> ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.bed

#
# Label annotation and simplify by flattening.
#
# Add 'ASD_Symbol=' prefix and
# remove duplicate values in one go using perl:
#
perl -F'\t' -lane \
  '$F[3] = join(",ASD_Symbol=", sort(grep(!$_{$_}++, split(",", $F[3])))); 
   print join "\t", @F; %_ = ();
  ' \
  ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.bed \
| awk -F "\t" 'BEGIN {OFS=FS} {print $1,$2,$3,"ASD_Symbol="$4}' \
> ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
mv ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed{,.tmp}
for char_class in {'0-9','X','Y','M'}; do \
  awk "\$1 ~ /^[${char_class}]/" \
  ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp \
  | sort -k1V -k2n \
  >> ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed; \
done
rm ${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp

#################################################################################
# 4.  Compare Ensembl and ASD BED files and merge into single CGD BED file
#################################################################################

cd ${CGD_DIR}/BED/
export ENSEMBL_BED_FILE=${CGD_DIR}/BED/From_Ensembl_API/${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
export     ASD_BED_FILE=${CGD_DIR}/BED/From_ASD/${CGD_FILE}.${TARGET}-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

#
# Count total number of regions:
#
wc -l \
${ENSEMBL_BED_FILE} \
${ASD_BED_FILE}
#
#  46752 /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_Ensembl_API/CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#  46507 /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_ASD/CGD_2014-10-08_updated.txt.all-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#

#
# Count total number of covered bases:
#
for bed_file in {${ENSEMBL_BED_FILE},${ASD_BED_FILE}}; do \
   awk '{total_bases+=$3-$2} END {printf "   %s", total_bases}' ${bed_file}; \
   echo " ${bed_file}";
done
#
#   18536401 /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_Ensembl_API/CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#   17309686 /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_ASD/CGD_2014-10-08_updated.txt.all-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#

#
# Compare Ensembl to ASD
#
# DO NOT use -sorted option.
# Our data is sorted, just not in the order intersectBed expects it in resulting in many many false positive differences...
#
ensembl_only_regions=$(intersectBed -v -a ${ENSEMBL_BED_FILE} -b ${ASD_BED_FILE} | wc -l | sed 's/ *//')
asd_only_regions=$(intersectBed -v -b ${ENSEMBL_BED_FILE} -a ${ASD_BED_FILE} | wc -l | sed 's/ *//')
ensembl_only_bases=$(subtractBed -a ${ENSEMBL_BED_FILE} -b ${ASD_BED_FILE} | awk '{total_bases+=$3-$2} END {printf "%s", total_bases}')
asd_only_bases=$(subtractBed -b ${ENSEMBL_BED_FILE} -a ${ASD_BED_FILE} | awk '{total_bases+=$3-$2} END {printf "%s", total_bases}')
echo "==============================================================
Counted ${ensembl_only_regions} regions only present in ${ENSEMBL_BED_FILE}
Counted ${asd_only_regions} regions only present in ${ASD_BED_FILE}
==============================================================
Counted ${ensembl_only_bases} bases only present in ${ENSEMBL_BED_FILE}
Counted ${asd_only_bases} bases only present in ${ASD_BED_FILE}
=============================================================="
#
# ==============================================================
# Counted 1801 regions only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_Ensembl_API/CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
# Counted 757 regions only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_ASD/CGD_2014-10-08_updated.txt.all-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
# ==============================================================
# Counted 1857716 bases only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_Ensembl_API/CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
# Counted 631001 bases only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/From_ASD/CGD_2014-10-08_updated.txt.all-targetsFromASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
# ==============================================================
#

#
# Create merged, sorted BED file.
#
export MERGED_CGD_BED_FILE_PREFIX="${CGD_FILE}.${TARGET}-targetsFromEnsembl${ENSEMBL_VERSION}AndASD"
cat ${ENSEMBL_BED_FILE} ${ASD_BED_FILE} \
| sort -k1,1V -k2,2n \
| mergeBed -nms \
  -i stdin \
> ${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed

#
# Sort annotation and remove duplicates
#
perl -F'\t' -lane'
  $F[3] = join(",", sort(grep(!$_{$_}++, split(",", $F[3])))); 
  print join "\t", @F; %_ = ();
  ' \
  ${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed \
> ${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
mv ${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed{,.tmp}
for char_class in {'0-9','X','Y','M'}; do \
  awk "\$1 ~ /^[${char_class}]/" \
  ${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp \
  | sort -k1V -k2n \
  >> ${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed; \
done
rm ${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp

#################################################################################
# 5.  Compare and merge into one BED file:
#     * CGD 2014-10-08
#     * Agilent SureSelect Inherited Disease V1 capturing kit (ELID S0684402), 
#       which was created as described in:
#       ${ROOT_BED_FILE_DIR}/Agilent_SureSelect_Inherited_Disease_V1_S0684402/ \
#       README_create_Agilent_SureSelect_Inherited_Disease_V1_S0684402_BED_file.txt
#################################################################################

#
# DO NOT use -sorted option.
# Our data is sorted, just not in the order intersectBed expects it in resulting in many many false positive differences...
#
export SSID_BED=${ROOT_BED_FILE_DIR}/Agilent_SureSelect_Inherited_Disease_V1_S0684402/gcc_version_indexed_for_human_g1k_v37/\
Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.labelled.bed
export MERGED_CGD_BED=${CGD_DIR}/BED/${MERGED_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

#
# Count total number of covered bases:
#
for bed_file in {${SSID_BED},${MERGED_CGD_BED}}; do \
   awk '{total_bases+=$3-$2} END {printf "   %s", total_bases}' ${bed_file}; \
   echo " ${bed_file}";
done
#
#    8.343.358 /gcc/sources/BED_Files//Agilent_SureSelect_Inherited_Disease_V1_S0684402/gcc_version_indexed_for_human_g1k_v37/Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.labelled.bed
#   19.167.402 /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75AndASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#

ssid_only_regions=$(intersectBed -v -a ${SSID_BED} -b ${MERGED_CGD_BED} | wc -l | sed 's/ *//')
cgd_only_regions=$(intersectBed  -v -b ${SSID_BED} -a ${MERGED_CGD_BED} | wc -l | sed 's/ *//')
ssid_only_bases=$(subtractBed       -a ${SSID_BED} -b ${MERGED_CGD_BED} | awk '{total_bases+=$3-$2} END {printf "%s", total_bases}')
cgd_only_bases=$(subtractBed        -b ${SSID_BED} -a ${MERGED_CGD_BED} | awk '{total_bases+=$3-$2} END {printf "%s", total_bases}')
echo "==============================================================
Counted ${ssid_only_regions} regions only present in ${SSID_BED}
Counted ${cgd_only_regions} regions only present in ${MERGED_CGD_BED}
==============================================================
Counted ${ssid_only_bases} bases only present in ${SSID_BED}
Counted ${cgd_only_bases} bases only present in ${MERGED_CGD_BED}
=============================================================="
#
# ==============================================================
# Counted  4.151 regions only present in /gcc/sources/BED_Files//Agilent_SureSelect_Inherited_Disease_V1_S0684402/gcc_version_indexed_for_human_g1k_v37/Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.labelled.bed
# Counted 12.490 regions only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75AndASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
# ==============================================================
# Counted    851.966 bases only present in /gcc/sources/BED_Files//Agilent_SureSelect_Inherited_Disease_V1_S0684402/gcc_version_indexed_for_human_g1k_v37/Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.labelled.bed
# Counted 11.676.010 bases only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75AndASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
==============================================================
#
# Unique bases in CGD     compared to size of SSID V1: 11.676.010 /  8.343.358 * 100 =~ 139.9 %
# Unique bases in SSID V1 compared to size of CGD:        851.966 / 19.167.402 * 100 =~   4.4 %
#

#
# Create merged Agilent SureSelect Inherited Disease + CGD BED file
#
export MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX="Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated"
cat ${SSID_BED} ${MERGED_CGD_BED} \
| sort -k1V -k2n \
| mergeBed -nms \
  -i stdin \
> ${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed

#
# Sort annotation and remove duplicates
#
perl -F'\t' -lane'
  $F[3] = join(",", sort(grep(!$_{$_}++, split(",", $F[3])))); 
  print join "\t", @F; %_ = ();
  ' \
  ${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed \
> ${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
mv ${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed{,.tmp}
for char_class in {'0-9','X','Y','M'}; do \
  awk "\$1 ~ /^[${char_class}]/" \
  ${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp \
  | sort -k1V -k2n \
  >> ${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed; \
done
rm ${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp

#
# Check size of merged BED file: covered bases shouls have increased compared to individual files.
#
for bed_file in {${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed,}; do \
   awk '{total_bases+=$3-$2} END {printf "   %s", total_bases}' ${bed_file}; \
   echo " ${bed_file}";
done
#
#   20.019.368 /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#

#################################################################################
# 6.  Compare and merge into one BED file:
#     * "Clinical Exome": CGD 2014-10-08 + Agilent SureSelect Inherited Disease V1 capturing kit (ELID S0684402).
#     * "Research Exome": Agilent SureSelect Human All Exon V5 capturing kit (ELID S04380110).
#################################################################################

#
# DO NOT use -sorted option.
# Our data is sorted, just not in the order intersectBed expects it in resulting in many many false positive differences...
#
export EXOME_PADDED_BED_FILE_PREFIX=Agilent_SureSelect_Human_All_Exon_V5_S04380110_Padded
export EXOME_EXTRAPADDED_BED=${ROOT_BED_FILE_DIR}/Agilent_SureSelect_Human_All_Exon_V5_S04380110/gcc_version_indexed_for_human_g1k_v37/${EXOME_PADDED_BED_FILE_PREFIX}.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
export MERGED_SSID_PLUS_CGD_BED=${CGD_DIR}/BED/${MERGED_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

research_exome_only_regions=$(intersectBed -v -a ${EXOME_EXTRAPADDED_BED} -b ${MERGED_SSID_PLUS_CGD_BED} | wc -l | sed 's/ *//')
clinical_exome_only_regions=$(intersectBed -v -b ${EXOME_EXTRAPADDED_BED} -a ${MERGED_SSID_PLUS_CGD_BED} | wc -l | sed 's/ *//')
research_exome_only_bases=$(subtractBed -a ${EXOME_EXTRAPADDED_BED} -b ${MERGED_SSID_PLUS_CGD_BED} | awk '{total_bases+=$3-$2} END {printf "%s", total_bases}')
clinical_exome_only_bases=$(subtractBed -b ${EXOME_EXTRAPADDED_BED} -a ${MERGED_SSID_PLUS_CGD_BED} | awk '{total_bases+=$3-$2} END {printf "%s", total_bases}')
echo "==============================================================
Counted ${research_exome_only_regions} regions only present in ${EXOME_EXTRAPADDED_BED}
Counted ${clinical_exome_only_regions} regions only present in ${MERGED_SSID_PLUS_CGD_BED}
==============================================================
Counted ${research_exome_only_bases} bases only present in ${EXOME_EXTRAPADDED_BED}
Counted ${clinical_exome_only_bases} bases only present in ${MERGED_SSID_PLUS_CGD_BED}
=============================================================="
#
# ==============================================================
# Counted 142.676 regions only present in /gcc/sources/BED_Files//Agilent_SureSelect_Human_All_Exon_V5_S04380110/gcc_version_indexed_for_human_g1k_v37/Agilent_SureSelect_Human_All_Exon_V5_S04380110_Padded.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
# Counted   6.011 regions only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
# ==============================================================
# Counted 78.270.215 bases only present in /gcc/sources/BED_Files//Agilent_SureSelect_Human_All_Exon_V5_S04380110/gcc_version_indexed_for_human_g1k_v37/Agilent_SureSelect_Human_All_Exon_V5_S04380110_Padded.sorted.stripped.flanked0-20bp.regionMerged.labelled.bed
# Counted  8.722.572 bases only present in /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
# ==============================================================
#
# Unique bases in Clinical Exome compared to size of Research Exome:  8.722.572 / 89.567.011 * 100 =~   9.7 %
# Unique bases in Research Exome compared to size of Clinical Exome: 78.270.215 / 20.019.368 * 100 =~ 391.0 %
#

#
# Create merged Agilent SureSelect All Exon V5 capturing kit + CGD BED file
#
export MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX="Agilent_SureSelect_Human_All_Exon_V5_S04380110_Padded_and_Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated"
cat ${EXOME_EXTRAPADDED_BED} ${MERGED_SSID_PLUS_CGD_BED} \
| sort -k1V -k2n \
| mergeBed -nms \
  -i stdin \
> ${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed

#
# Sort annotation and remove duplicates
#
perl -F'\t' -lane'
  $F[3] = join(",", sort(grep(!$_{$_}++, split(",", $F[3])))); 
  print join "\t", @F; %_ = ();
  ' \
  ${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.bed \
> ${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed

#
# Sort for human_g1k_v37.fa reference genome used in GCC pipelines.
#
mv ${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed{,.tmp}
for char_class in {'0-9','X','Y','M'}; do \
  awk "\$1 ~ /^[${char_class}]/" \
  ${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp \
  | sort -k1V -k2n \
  >> ${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed; \
done
rm ${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed.tmp

#
# Check size of merged BED file: covered bases shouls have increased compared to individual files.
#
for bed_file in {${CGD_DIR}/BED/${MERGED_EXOME_PLUS_SSID_PLUS_CGD_BED_FILE_PREFIX}.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed,}; do \
   awk '{total_bases+=$3-$2} END {printf "   %s", total_bases}' ${bed_file}; \
   echo " ${bed_file}";
done
#
#   98.289.583 /gcc/sources/BED_Files//ClinicalGenomicsDatabase_2014-10-08//BED/Agilent_SureSelect_Human_All_Exon_V5_S04380110_Padded_and_Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#

#################################################################################
# Final BED files used in 5GPM:
#################################################################################
# 1. Clinical Genomics Database (CGD) only:
#    /gcc/sources/BED_Files/ClinicalGenomicsDatabase_2014-10-08/BED/
#    CGD_2014-10-08_updated.txt.all-targetsFromEnsembl75AndASD.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#
# 2. Agilent Sure Select Inherited Disease (SSID) only
#    /gcc/sources/BED_Files/Agilent_SureSelect_Inherited_Disease_V1_S0684402/gcc_version_indexed_for_human_g1k_v37/
#    Agilent_SureSelect_Human_Inherited_Disease_V1_S0684402_UMCG.sorted.stripped.flanked20bp.regionMerged.labelled.bed
#
# 3. Clinical Exome: CGD + SSID
#    /gcc/sources/BED_Files/ClinicalGenomicsDatabase_2014-10-08/BED/
#    Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#
# 4. Clinical + Research Exome: CGD + SSID + Agilent Sure Select Human All Exon:
#    /gcc/sources/BED_Files/ClinicalGenomicsDatabase_2014-10-08/BED/
#    Agilent_SureSelect_Human_All_Exon_V5_S04380110_Padded_and_Agilent_SureSelect_Inherited_Disease_V1_S0684402_UMCG_and_CGD_2014-03-12_updated.sorted.flanked20bp.regionMerged.flattenedAnnotation.bed
#################################################################################
