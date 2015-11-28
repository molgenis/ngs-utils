#!/usr/bin/env bash

#
# Bash sanity
#
set -e  # Exit if any subcommand or pipeline returns a non-zero status.
set -u  # Exit if any uninitialised variable is used.

usage() {
    echo '##################################################################################################'
    echo ' This tool compares 2 VCF files and shows only the different rows based on the following columns:'
    echo '   CHROM, POS, REF, ALT, FILTER, INFO:DP and SampleID:GT'
    echo ' The read depth (INFO:DP) can be used to select only variants with a read depth >= threshold.'
    echo ' This tool will work well if the amount of differences is small, '
    echo ' but may become slow when the number differences is large.'
    echo
    echo ' Usage:'
    echo '   bash vcf-compare.sh -1 <input_vcf_file> -2 <input_vcf_file> -d <minimal_depth> -o <output_folder>'
    echo
    echo ' Example:'
    echo '   bash vcf-compare.sh -1 old_analysis.snps.final.vcf \'
    echo '                       -2 new_analysis.snps.final.vcf \'
    echo '                       -d 20 \'
    echo '                       -o ./results/'
    echo '##################################################################################################'
}

declare VCF1=""
declare VCF2=""
declare DEPTH="0"
declare OUT=""

#
# Take the parameters given on the commandline and assign them to variables.
#
while getopts ":1:2:d:o:h" option; do
    case "${option}" in
        1)  VCF1=${OPTARG};;
        2)  VCF2=${OPTARG};;
        d) DEPTH=${OPTARG};;
        o)   OUT=${OPTARG};;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "ERROR: Invalid option -${OPTARG}. Try \"$(basename $0) -h\" for help."
            exit 1
            ;;
        :)
            echo "ERROR: Option -${OPTARG} requires an argument. Try \"$(basename $0) -h\" for help."
            ;;
    esac
done

#
# Check if all parameters are set.
#
if [[ ${VCF1} && ${VCF2} && ${OUT} ]]; then
    echo
    echo "Comparing ${VCF1} to ${VCF2}..."
    echo
else
    usage
    echo
    echo "ERROR: missing required argument. Try \"$(basename $0) -h\" for help."
    echo
    exit 1
fi

#
# Create tmp folders.
#
SCRATCH="${OUT}/TMP/"
rm -Rf   ${SCRATCH}
mkdir -p ${SCRATCH}

#
# select filenames from given path and check if they are unique.
#
vcf01="$(basename ${VCF1})"
vcf02="$(basename ${VCF2})"
if [ ${vcf01} == ${vcf02} ]; then
	usage
	echo
	echo "ERROR: ${vcf01} is equal to ${vcf02}."
	echo "       Make sure the filenames are unique!"
	echo
	exit 1
fi

#
# We need VCFtools.
#
echo
echo "Locating vcftools..."
module load vcftools
module list
echo



#
# Extract depth (DP) from INFO field dropping any other INFO subfields.
#
vcftools --recode --recode-INFO DP --vcf ${VCF1} --out "${SCRATCH}/${vcf01}.stripped"
vcftools --recode --recode-INFO DP --vcf ${VCF2} --out "${SCRATCH}/${vcf02}.stripped"
bgzip "${SCRATCH}/${vcf01}.stripped.recode.vcf"
bgzip "${SCRATCH}/${vcf02}.stripped.recode.vcf"
tabix -p vcf "${SCRATCH}/${vcf01}.stripped.recode.vcf.gz"
tabix -p vcf "${SCRATCH}/${vcf02}.stripped.recode.vcf.gz"

#
# Select only the following columns: CHROM (1), POS (2), REF(4), ALT(5), FILTER(7), INFO:DP and SampleID:GT (10). 
#
vcf-query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DP[\t%GT]\n' "${SCRATCH}/${vcf01}.stripped.recode.vcf.gz" > "${SCRATCH}/${vcf01}.stripped.txt"
vcf-query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DP[\t%GT]\n' "${SCRATCH}/${vcf02}.stripped.recode.vcf.gz" > "${SCRATCH}/${vcf02}.stripped.txt"
#vcf-query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DP[\t%SAMPLE=%GT]\n' "${OUT}/TMP01/${vcf01}.stripped.recode.vcf" > "${OUT}/TMP01/${vcf01}.stripped.txt"
#vcf-query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DP[\t%SAMPLE=%GT]\n' "${OUT}/TMP02/${vcf02}.stripped.recode.vcf" > "${OUT}/TMP02/${vcf02}.stripped.txt"

if [ ${DEPTH} > 0 ]; then
    #
    # Remove variants with DP <= $DEPTH. 
    #
    awk -F "\t" 'BEGIN { OFS=FS } $6 > '"${DEPTH}"'-1 { print $1,$2,$3,$4,$5,$6,$7 }' "${SCRATCH}/${vcf01}.stripped.txt" > "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt"
    awk -F "\t" 'BEGIN { OFS=FS } $6 > '"${DEPTH}"'-1 { print $1,$2,$3,$4,$5,$6,$7 }' "${SCRATCH}/${vcf02}.stripped.txt" > "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt"
else
    awk -F "\t" 'BEGIN { OFS=FS } { print $1,$2,$3,$4,$5,$6,$7 }' "${SCRATCH}/${vcf01}.stripped.txt" > "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt"
    awk -F "\t" 'BEGIN { OFS=FS } { print $1,$2,$3,$4,$5,$6,$7 }' "${SCRATCH}/${vcf02}.stripped.txt" > "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt"
fi

#
# Check if we have any variants left.
#
if [[ ! -s "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt" ]]; then
	echo
	echo "WARN: ${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt is empty."
	echo "      There were no variants left after filtering for a minimal read depth of ${DEPTH}."
	echo
fi
if [[ ! -s "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt" ]]; then
	echo
	echo "WARN: ${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt is empty."
	echo "      There were no variants left after filtering for a minimal read depth of ${DEPTH}."
	echo
fi
if [[ ! -s "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt" && ! -s "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt" ]]; then
	echo
	echo "ERROR: There were no variants left after filtering for a minimal read depth of ${DEPTH} in either VCF file."
	echo
	exit 1
fi 

#
# Reformat and add header: create small output for diff.
#
vcf01_sample_IDs=$(vcf-query -l "${SCRATCH}/${vcf01}.stripped.recode.vcf.gz")
vcf02_sample_IDs=$(vcf-query -l "${SCRATCH}/${vcf02}.stripped.recode.vcf.gz")
printf "%-6s  %10s  %-10s  %-10s  %-s\n" "#CHROM" "POS" "REF" "ALT" "${vcf01_sample_IDs}:GT" \
  > "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt.small.reformatted"
awk -F "\t" '{printf "%-6s  %10s  %-10s  %-10s  %-s\n", $1,$2,$3,$4,$7 }' \
    "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt" \
 >> "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt.small.reformatted"
#
printf "%-6s  %10s  %-10s  %-10s  %-s\n" "#CHROM" "POS" "REF" "ALT" "${vcf02_sample_IDs}:GT" \
  > "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt.small.reformatted"
awk -F "\t" '{printf "%-6s  %10s  %-10s  %-10s  %-s\n", $1,$2,$3,$4,$7 }' \
    "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt" \
 >> "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt.small.reformatted"

#
# Intermediate cleanup.
#
mv ${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt.small{.reformatted,}
mv ${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt.small{.reformatted,}

#
# Intermediate comparison of the subset of VCF data for the two files using diff.
#  * List only lines NOT starting with "#". Hence skip all meta-data and header lines.
#  * List all lines starting with a [<>] and write them to a new file.
#    These will then be used to search the corresponding complete/original files for the same variants, 
#    so we have more annotation/data available to compare the variants.
#
vcf_diff_result="${SCRATCH}/diff_${vcf01}_vs_${vcf02}.txt.small"
diff -w -d \
   "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt.small" \
   "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt.small" \
 | grep -v '^#' \
 | grep -P '^[<>]' \
 | awk 'BEGIN {OFS="\t"} {print $2,$3}' \
 | sort -u -k1,1n -k2,2n \
 > "${vcf_diff_result}"

#
# Count the total amount of sites taken into account: same in both VCFs + only present in VCF 1 + only present in VCF2.
#
number_of_sites_total=$(cat \
    "${SCRATCH}/${vcf01}.stripped.minDP${DEPTH}.txt.small" \
    "${SCRATCH}/${vcf02}.stripped.minDP${DEPTH}.txt.small" \
 | grep -v '^#' \
 | awk '{print $1,$2}' \
 | sort -u \
 | wc -l)

#
# Get the detailed info for these variants from the more complete files.
#
#IFS="\n" read -a variant_seq_coordinates <<< $(cat "${vcf_diff_result}")
while read -r chr pos; do
    #echo "DEBUG: Searching for variant @ coordinate: ${chr}:${pos}."
    set +e
#    grep_match=$(grep -P "^${chr}\t${pos}\t" "${SCRATCH}/${vcf01}.stripped.txt")
#    grep_status=$?
#    if [ ${grep_status} -eq 0 ]; then
#        echo "$grep_match" >> "${SCRATCH}/${vcf01}.stripped.txt.for.report"
#    elif [ ${grep_status} -eq 1 ]; then
#        # No matches were found: write that explicitly to result to prevent spurious partial diff matches.
#        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${chr}" "${pos}" "-" "-" "." "." "./." >> "${SCRATCH}/${vcf01}.stripped.txt.for.report"
#    else
#        echo
#        echo "ERROR: Cannot grep \"^${chr}\\t${pos}\\t\" in ${SCRATCH}/${vcf01}.stripped.txt."
#        echo "       Exit value is: $?"
#        echo
#        exit 1
#    fi
#    grep_match=$(grep -P "^${chr}\t${pos}\t" "${SCRATCH}/${vcf02}.stripped.txt")
#    grep_status=$?
#    if [ ${grep_status} -eq 0 ]; then
#        echo "$grep_match" >> "${SCRATCH}/${vcf02}.stripped.txt.for.report"
#    elif [ ${grep_status} -eq 1 ]; then
#        # No matches were found: write that explicitly to result to prevent spurious partial diff matches.
#        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${chr}" "${pos}" "-" "-" "." "." "./." >> "${SCRATCH}/${vcf02}.stripped.txt.for.report"
#    else
#        echo
#        echo "ERROR: Cannot grep \"^${chr}\\t${pos}\\t\" in ${SCRATCH}/${vcf02}.stripped.txt."
#        echo "       Exit value is: $?"
#        echo
#        exit 1
#    fi
    query_match=$(vcf-query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DP[\t%GT]\n' -r ${chr}:${pos}-${pos} "${SCRATCH}/${vcf01}.stripped.recode.vcf.gz")
    if [ -n "${query_match}" ]; then
        echo "$query_match" >> "${SCRATCH}/${vcf01}.stripped.txt.for.report"
    else
        # No matches were found: write that explicitly to result to prevent spurious partial diff matches.
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${chr}" "${pos}" "-" "-" "." "." "./." >> "${SCRATCH}/${vcf01}.stripped.txt.for.report"
    fi
    query_match=$(vcf-query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/DP[\t%GT]\n' -r ${chr}:${pos}-${pos} "${SCRATCH}/${vcf02}.stripped.recode.vcf.gz")
    if [ -n "${query_match}" ]; then
        echo "$query_match" >> "${SCRATCH}/${vcf02}.stripped.txt.for.report"
    else
        # No matches were found: write that explicitly to result to prevent spurious partial diff matches.
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${chr}" "${pos}" "-" "-" "." "." "./." >> "${SCRATCH}/${vcf02}.stripped.txt.for.report"
    fi
    set -e
done < "${vcf_diff_result}"

#
# Create report with diff's side-by-side comparison feature.
#
vcf_diff_report="${SCRATCH}/diff_${vcf01}_vs_${vcf02}.txt"
diff -y -w -W 150 --horizon-lines 1 \
   "${SCRATCH}/${vcf01}.stripped.txt.for.report" \
   "${SCRATCH}/${vcf02}.stripped.txt.for.report" \
 | grep -P '[\|\<\>]' \
 > "${vcf_diff_report}"
number_of_sites_discordant=$(wc -l ${vcf_diff_report} | awk '{print $1}')

#
# Reformat and add header: full output for extra annotation to inspect diff results.
#
vcf01_sample_IDs_length=$(echo "${#vcf01_sample_IDs} + 3" | bc)
vcf02_sample_IDs_length=$(echo "${#vcf02_sample_IDs} + 3" | bc)
format="%-6s  %10s  %-10s  %-10s  %-8s  %-7s  %-${vcf01_sample_IDs_length}s  |  %-6s  %10s  %-10s  %-10s  %-8s  %-7s  %-${vcf02_sample_IDs_length}s\n"
printf "${format}" \
       "#CHROM" "POS" "REF" "ALT" "FILTER" "INFO:DP" "${vcf01_sample_IDs}:GT" \
       "#CHROM" "POS" "REF" "ALT" "FILTER" "INFO:DP" "${vcf02_sample_IDs}:GT" \
    > "${vcf_diff_report}.reformatted"
awk -v format="%-6s  %10s  %-10s  %-10s  %-8s  %-7s  %-${vcf01_sample_IDs_length}s  %1s  %-6s  %10s  %-10s  %-10s  %-8s  %-7s  %-${vcf02_sample_IDs_length}s\n" \
    '{printf format, $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15 }' \
    "${vcf_diff_report}" \
 >> "${vcf_diff_report}.reformatted"

#
# Intermediate cleanup.
#
final_result="${OUT}/diff_${vcf01}_vs_${vcf02}.txt"
mv "${vcf_diff_report}.reformatted" "${final_result}"

#
# Calculate discordance.
#
discordance=$(echo "scale=6; ${number_of_sites_discordant}/${number_of_sites_total}*100" | bc)
discordance=$(printf "%.2f" ${discordance})

#
# Display result and stats.
#
line_length=$(echo "132 + ${vcf01_sample_IDs_length} + ${vcf02_sample_IDs_length}" | bc)
echo
head -c ${line_length} </dev/zero | tr '\0' '='; echo
cat "${final_result}"
head -c ${line_length} </dev/zero | tr '\0' '='; echo
echo
echo "INFO: For ${number_of_sites_total} sites where the variant call had a depth of at least ${DEPTH} in at least one of the VCF files:"
echo "INFO: Found ${number_of_sites_discordant}/${number_of_sites_total} sites (${discordance}%) with different variant calls (ignoring VCF FILTER values)."
echo "INFO: Finished comparison."
echo "INFO: The result displayed above is in ${final_result}."

exit 0
