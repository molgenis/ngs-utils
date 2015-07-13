# Clinival Pathogenicity Comparison

The scripts presented in this can be used to create a single merged dataset using different population datatsets.
The resulting dataset can be annotated with different kinds of extra information from external sources, for example ClinVar data or CADD scores.

The procedure below tells you how to use these scripts and create your very own merged and annotated population dataset in VCF format.

# dependecies
The dependencies listed below are need to run the all the scripts, each script has a usage shown when executed without arguments.  
 - Ensembl API : www.ensembl.org/info/docs/api/core/index.html#api
 - Vcftools : http://vcftools.sourceforge.net/downloads.html
 - Tabix : http://www.htslib.org/doc/tabix.html
 - Bgzip : http://www.htslib.org/doc/tabix.html


## 1. Datasets
To start off using these scripts it is needed to get the following datasets.
 - GoNL : https://molgenis26.target.rug.nl/downloads/gonl_public/variants/release5_with_GTC/release5_noContam_noChildren_with_AN_AC_GTC_stripped.tgz + X from release 4
 - 1000 genomes : http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
 - ExAC : ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/
 - ClinVar : ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz + ../vcf.gz/tbi
 - CADD : http://cadd.gs.washington.edu/download version 1.2
 - FitCon : http://cadd.gs.washington.edu/download All possible SNVs of GRCh37/hg19 incl. all annotations
 - DANN : https://cbcl.ics.uci.edu/public_data/DANN/data/

## 2. Merging
 By using any merging tool it is possible to create one big file to be used as a template for annotating.
 In the merge folder a script is provided for merging different datasets based on the VCF standard (http://samtools.github.io/hts-specs/VCFv4.3.pdf), however this script will strip away the INFO collumn to ensure a less verbose input file.

## 3. Filtering
 By using the provided CanonicalSeqRetriever script it is possible to get all exonic regions from genome build grch37 (http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/). 
 The script uses the Ensembl API becouse of this it cannot function without the API being installed. The resulting bed file (canonical_exon_sequences_filtered_headered_sorted.bed) is also located in the Bedfiles folder and has already been sorted and headered and ready for use.
 
 With the use of vcftools you can use this bed file (or any bed file) to extract the locations of interest.
 There is a bed file filter script generator available to create bash scripts to conduct this filtering.
 Besides the initial filtering it also compresses the file using bgzip and creates tabix index files for the vcf files.
 
## 4. Annotation
 With the use of the molgenis annotators (https://github.com/molgenis/molgenis) it is possible to add the information of other resources to the desired dataset.
 The commandline version only supports VCF files which is provided in annotator folder. 
 


 
