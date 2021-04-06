#!/bin/bash

#
##
### Hardcoded variables.
##
#
FASTA_FILE_EXTENSION='fa'

#
# Source URLs.
#
#GRCh38_FULL_AS='ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz'
#hs37d5='ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'

declare -a REFGENOMES=('GRCh38_HMF_AS1')

declare -A REFGENOME_DETAILS=(
    #
    # GRCh38 flavors
    #
    ['GRCh38_HMF_AS1,description']='Human Genome Analysis Set for the Hartwig Medical Foundation (release 1). Primary assembly unit based on GRCh38.'
    ['GRCh38_HMF_AS1,ingredients']='GRCh38_NoALT_AS hs38d1 PhiX174'
    ['GRCh38_HMF_AS1,md5']='6e7f8a01b7c7cb742b323c9df28a40fa'
)

declare -A ALL_INGREDIENTS=(
    ['GRCh38_NoALT_AS,description']='GRCh38 analysis set including chromosomes, EBV viral genome, unplaced contigs and unlocalized contigs. (Excluding ALT and PATCH contigs.)'
    #['GRCh38_NoALT_AS,url']='ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    #['GRCh38_NoALT_AS,url']='ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    #['GRCh38_NoALT_AS,url']='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    ['GRCh38_NoALT_AS,url']='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz'
    ['GRCh38_NoALT_AS,md5']='beadb3e9633f6a39bbf50e3ba995468c'
    ['hs38d1,description']='Human decoy sequences for GRCh38.'
    #['hs38d1,url']='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz'
    ['hs38d1,url']='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/786/075/GCA_000786075.2_hs38d1/GCA_000786075.2_hs38d1_genomic.fna.gz'
    ['hs38d1,md5']='b7e2a3f3193a2d82d0f4bea7da4de398'
    ['PhiX174,description']='PhiX174 genome (used as decoy for spike-in controls).'
    #['PhiX174,url']='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz'
    ['PhiX174,url']='ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz'
    ['PhiX174,md5']='91d568118a1e14f1105a79e71928dc31'
)

declare -A REFGENOME_READMES=()
REFGENOME_READMES['GRCh38_HMF_AS1']="
==============================================================================
 README for GRCh38_HMF_AS1
==============================================================================

This is a Human Genome Analysis Set for the Hartwig Medical Foundation (release 1).
The primary assembly unit based is on GRCh38 (GCA_000001405.15), 
but the 'analysis set' has a few modifications compared to the orginal reference genome 
to optimise for use by various Next Generation Sequence analysis pipelines.

==============================================================================
 The ingredients:
==============================================================================

A. ${ALL_INGREDIENTS['GRCh38_NoALT_AS,url']}

   A gzipped FastA file containing the following sequences:
1. Autosomal and sex chromosomes from the GRCh38 Primary Assembly unit.
   * The two PAR regions on chrY have been hard-masked with Ns. 
     The chromosome Y sequence provided therefore has the same 
     coordinates as the original sequence but it is not identical.
   * Duplicate copies of centromeric arrays and WGS 
     on chromosomes 5, 14, 19, 21 & 22 have been hard-masked with Ns. 
2. Mitochondrial genome from the GRCh38 non-nuclear assembly unit.
   This is the revised Cambridge Reference Sequence a.k.a rCRS (NC_012920.1)
3. Unlocalized scaffolds from the GRCh38 Primary Assembly unit.
4. Unplaced scaffolds from the GRCh38 Primary Assembly unit.
5. Epstein-Barr Virus (EBV) sequence.
   The EBV sequence is not part of the human genome assembly, 
   but is included as a decoy for alignment of reads, 
   that are often present in sequencing samples either 
   as result of a 'natural infection' or because EBV was used to immortalise cell lines.

B. ${ALL_INGREDIENTS['hs38d1,url']}

   A gzipped FastA file containing the following sequences:
6. Human decoy sequences from hs38d1 (GCA_000786075.2)

C. ${ALL_INGREDIENTS['PhiX174,url']}

   A gzipped FastA file containing the following sequences:
7. The phage PhiX174 genome (NC_001422.1)
   The PhiX174 sequence is not part of the human genome assembly
   but is included as a decoy for alignment of reads, 
   that are often present in sequencing samples, 
   because they were added as in vitro or in silico spike-in controls.
   
==============================================================================
 FastA sequence meta-data lines
==============================================================================

All FastA meta-data lines were patched to remove unstable IDs, create uniformity 
and enable easy reporting of variants in HGVS compliant format.
All meta-data after the FastA sequence identifier (averything after the > until the first white space)
is in the format of double space separated key:value pairs.
Some examples:

>NC_000024.10(Y)  AC:CM000686.2  LN:57227415  rl:Chromosome  M5:ce3e31103314a704255f3cd90369ecce  AS:GRCh38  hm:10001-2781479,56887903-57217415
>NC_012920.1(MT)  AC:J01415.2  LN:16569  rl:Mitochondrion  M5:c68f52674c9fb33aef52dcf399755519  AS:GRCh38  tp:circular
>NC_007605.1(EBV)  AC:AJ507799.2  LN:171823  rl:decoy  M5:6743bd63b3ff2b5b8985d8933c53290a  SP:Human_herpesvirus_4  tp:circular
>KN707606.1(UN_DECOY_00416)  AC:KN707606.1  LN:2200  rl:decoy  M5:20c768ac79ca38077e5012ee0e5f8333  AS:hs38d1
>NC_001422.1(PhiX174)  AC:J02482.1  LN:5386  rl:decoy  M5:3332ed720ac7eaa9b3655c06f6b9e196  SP:Enterobacteria_phage_phiX174_sensu_lato  tp:circular

Tag  Description
>    Versioned accession number preferably from RefSeq followed by a sequence name between round brackets.
     This allows for easy creation of HGVS complient variants like this NC_012920.1(MT):g.8860A>G.
     Note that although MT:g.8860A>G is also HGVS compliant reporting variants this way is useless 
     as it is unclear which version of which sequence and therefore which coordinate system was used.
AC:  Accession number of primary sequence database = INSD = ENA/GenBank/DDBJ.
LN:  Sequence length
rl:  Role of the sequence in this reference genome. One of:
       Chromosome
       Mitochondrion
       decoy
       unlocalized (Chromosome known, but exact position unknown.)
       unplaced    (Chromosome unknown)
M5  MD5 checksum of the sequence without the header / meta-data line and as one large line without newline characters.
AS  Assembly name.
hm  Hard Masked regions. Either a single span, two spans separated by a comma or 'multiple'.
tp  Topology. Either 'circular' for MT, EBV and PhiX or not specified for linear sequences.
"
#
##
### Functions.
##
#
function showHelp() {
  #
  # Display commandline help on STDOUT.
  #
  cat <<EOH

Usage:

   $(basename $0) -l|c [-k] -r genome -o /path/to/results/

Details:

   -r genome : Compose reference genome.
   -l        : List ingredients for the genome specified with -r.
   -c        : Compose the genome specified with -r.
   -k        : Keep the seperate ingredients (optional).
   -o dir    : Output directory. Only required for compiling genomes.

Reference genomes known to this script: 

EOH

for GENOME in "${REFGENOMES[@]}"; do
  printf "   %-16s : %s\n" "${GENOME}" "${REFGENOME_DETAILS["${GENOME},description"]}"
done

echo

}

function elementExistsIn () {
  for _ELEMENT in "${@:2}"; do
    [[ "${_ELEMENT}" == "${1}" ]] && return 0;
  done;
  return 1
}

function listIngredients() {
  local    _GENOME="${1}"
  local -a _MY_INGREDIENTS=("${REFGENOME_DETAILS["${_GENOME},ingredients"]}")
  local    _i=''
  echo "INFO: Reference Genome ${_GENOME} contains the following ingredients:"
  for _i in ${_MY_INGREDIENTS[@]}; do
    echo "INFO:  * ${_i}"
    local _MY_URL="${ALL_INGREDIENTS["${_i},url"]}"
    local _MY_DES="${ALL_INGREDIENTS["${_i},description"]}"
    echo "INFO:      Source:      ${_MY_URL}"
    echo "INFO:      Description: ${_MY_DES}"
  done
}

function compileReferenceGenome() {
  local    _GENOME="${1}"
  local    _RESULT_DIR="${2}"
  local -a _MY_INGREDIENTS=("${REFGENOME_DETAILS["${_GENOME},ingredients"]}")
  local    _i=''
  #
  # Check if FastA exists to prevent accidentally overwriting a reference genome.
  #
  if [ -f "${_RESULT_DIR}/${_GENOME}.${FASTA_FILE_EXTENSION}" ]; then
    echo "FATAL: FastA file for reference genome ${_GENOME} already exists and I will refuse to overwrite it."
    exit 1
  fi
  #
  # Try to create result dir if it is not already present and check if we can write there if it is present.
  #
  if [ ! -d "${_RESULT_DIR}" ]; then
    mkdir -p "${_RESULT_DIR}"
  fi
  if [[ ! -r "${_RESULT_DIR}" || ! -w "${_RESULT_DIR}" || ! -x "${_RESULT_DIR}" ]]; then
    echo "FATAL: You do not have enough permissions to store the result in dir: ${_RESULT_DIR}"
    exit 1
  fi
  #
  # Loop over list of ingredients and fetch, patch, verify checksum & append.
  #
  for _i in ${_MY_INGREDIENTS[@]}; do
    local _MY_URL="${ALL_INGREDIENTS["${_i},url"]}"
    local _MY_MD5="${ALL_INGREDIENTS["${_i},md5"]}"
    echo "INFO: Processing ingredient ${_i}..."
    echo "INFO:     Fetching ingredient from:"
    echo "          ${_MY_URL}..."
    (wget -O- ${_MY_URL} | gzip -dc) > "${_RESULT_DIR}/${_i}.${FASTA_FILE_EXTENSION}"
    if [ -r "${PATCHES_DIR}/${_i}.patch" ]; then
      echo "INFO:     Applying patch ${PATCHES_DIR}/${_i}.patch to ${_i}.${FASTA_FILE_EXTENSION}"
      patch "${_RESULT_DIR}/${_i}.${FASTA_FILE_EXTENSION}" "${PATCHES_DIR}/${_i}.patch"
    fi
    printf '%s  %s\n' "${_MY_MD5}" "${_i}.${FASTA_FILE_EXTENSION}" > "${_RESULT_DIR}/${_i}.md5"
    echo -n "INFO:     Verifying MD5 checksum for "
    (cd "${_RESULT_DIR}" && md5sum -c "${_i}.md5")
    echo "INFO:     Appending ${_i}.${FASTA_FILE_EXTENSION} to ${_GENOME}.${FASTA_FILE_EXTENSION}"
    cat "${_RESULT_DIR}/${_i}.${FASTA_FILE_EXTENSION}" >> "${_RESULT_DIR}/${_GENOME}.${FASTA_FILE_EXTENSION}"
    if [ ${KEEP_INGREDIENTS} -eq 0 ]; then
      rm "${_RESULT_DIR}/${_i}.${FASTA_FILE_EXTENSION}"  "${_RESULT_DIR}/${_i}.md5"
    fi
  done
  #
  # Verify checksum and create README for this reference genome.
  #
  printf '%s  %s\n' "${REFGENOME_DETAILS["${_GENOME},md5"]}" "${_GENOME}.${FASTA_FILE_EXTENSION}" > "${_RESULT_DIR}/${_GENOME}.md5"
  echo -n "INFO: Verifying MD5 checksum for "
  (cd "${_RESULT_DIR}" && md5sum -c "${_GENOME}.md5")
  echo -n "INFO: Creating ${_GENOME}.README... "
  cat > "${_RESULT_DIR}/${_GENOME}.README"  <<EOR
${REFGENOME_READMES[${_GENOME}]}
EOR
  echo 'done.'
  #
  # Signal success.
  #
  echo "INFO: Successfully compiled ${GENOME} reference genome in ${_RESULT_DIR}/."
  #
  # Provide advise on how to index reference for use by various NGS analysis tools.
  #
  echo 'INFO: ======================================================================================'
  echo 'INFO:  In order to use your new Reference Genome with common NGS analysis pipelines '
  echo 'INFO:  you may have create indices and other derived meta-data using the code listed below:'
  echo 'INFO: ======================================================================================'
  cat <<EOA
       #
       # Load latest versions of samtools, Picard and BWA.
       #
       module load SAMtools
       module load picard
       module load BWA
       module list
       #
       # Use samtools to index FastA file.
       #
       samtools faidx "${_RESULT_DIR}/${_GENOME}.${FASTA_FILE_EXTENSION}"
       #
       # Use Picard to create dict file.
       #
       java -jar \${EBROOTPICARD}/picard.jar CreateSequenceDictionary \\
            REFERENCE="${_RESULT_DIR}/${_GENOME}.${FASTA_FILE_EXTENSION}" \\
               OUTPUT="${_RESULT_DIR}/${_GENOME}.dict"
       #
       # Index reference genome for alignment with BWA.
       #
       bwa index -a bwtsw "${_RESULT_DIR}/${_GENOME}.${FASTA_FILE_EXTENSION}"
EOA
  echo 'INFO: ======================================================================================'
  echo "INFO: I'm signing off; Happy NGS analysing! $(basename $0) out."

}

#
##
### Bash sanity and initialisation.
##
#

#
# Bash sanity.
#
set -u
set -e

#
# Get path to directory where this script is located.
#
SCRIPT_DIR=$(cd -P "$( dirname "$0" )" && pwd)
SCRIPT_NAME=$(basename "$0" .bash)
PATCHES_DIR="${SCRIPT_DIR}/${SCRIPT_NAME}.patches"
START_TS=$(date "+%Y-%m-%d-T%H%M")

#
##
### Process commandline arguments.
##
#

#
# Get commandline arguments.
#
KEEP_INGREDIENTS=0 # Default is to delete ingredients and keep only the final product (Reference Genome in FastA format + a README).
while getopts ":hlckr:o:" opt; do
  case $opt in
    h)
      showHelp
      ;;
    c)
      ACTION='compose'
      ;;
    l)
      ACTION='list'
      ;;
    r)
      GENOME="${OPTARG}"
      ;;
    k)
      KEEP_INGREDIENTS=1
      ;;
    o)
      OUTPUT_DIR="${OPTARG}"
      ;;
    \?)
      echo "FATAL: Invalid option -${OPTARG}. Try \"$(basename $0) -h\" for help."
      exit 1
      ;;
    :)
      echo "FATAL: Option -${OPTARG} requires an argument. Try \"$(basename $0) -h\" for help."
      exit 1
      ;;
  esac
done

#
# Make sure there are no extra arguments we did not expect nor need.
#
shift $(($OPTIND - 1))
if [ ! -z ${1:-} ]; then
  echo "FATAL: Invalid argument \"$1\". Try \"$(basename $0) -h\" for help."
  exit 1
fi

#
##
### Main.
##
#

if [ -z "${GENOME:-}" ]; then
  #
  # No genome specified.
  #
  showHelp
  exit 1
elif elementExistsIn "${GENOME}" "${REFGENOMES[@]}"; then
  if [ "${ACTION:-}" == 'list' ]; then
    listIngredients "${GENOME}"
  elif [ "${ACTION:-}" == 'compose' ]; then
    if [ -z "${OUTPUT_DIR:-}" ]; then
      echo "FATAL: No output dir specified. Must tell me where to store the results..."
      exit 1
    fi
    compileReferenceGenome "${GENOME}" "${OUTPUT_DIR}"
  else
    echo "FATAL: No action specified. Must specify either -l for listing ingredients or -c for compiling the specified reference genome."
    exit 1
  fi
else
  echo "FATAL: Unknown reference genome specified: ${GENOME}."
  echo "INFO: Execute me without arguments for commandline help including a list of valid reference genomes."
  exit 1
fi

#[ ! -f $1.fa.bwt ] && echo -e "\nPlease run 'bwa index $1.fa'...\n"
