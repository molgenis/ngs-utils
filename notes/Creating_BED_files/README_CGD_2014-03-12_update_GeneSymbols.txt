#
# ======================================================
# How to update gene symbols
# referenced in the Clinical Genomics Database (CGD)
# to official gene symbols as curated by the HGNC
#
# Last update: 2014 March 19
# ======================================================
#

CURRENT_DATE=$(date "+%Y-%m-%d")

#
# Fetch new CGD from research.nhgri.nih.gov/CGD
#
curl http://research.nhgri.nih.gov/CGD/download/txt/CGD.txt.gz > CGD_${CURRENT_DATE}.txt.gz
gunzip CGD_${CURRENT_DATE}.txt.gz

#
# Create genesymbol2HGNC.txt by uploading gene symbols from CGD to HGNC gene symbol checker
#
# Go to http://www.genenames.org/cgi-bin/symbol_checker
# Paste the gene symbols from CGD_${CURRENT_DATE}.txt
# Save as genelist2HGNC_${CURRENT_DATE}.txt 
#
# Fetch "custom download" with Entrez Gene IDs
#
# Go to http://www.genenames.org/cgi-bin/download
# Select columns
#   HGNC ID
#   Synonyms
#   Approved Symbol
#   Name Synonyms
#   Entrez Gene ID (Curated by the HGNC)
#   Approved Name
#   Previous Symbols
#   Chromosome
#   Status
#   Previous Names
#   Entrez Gene ID (supplied by NCBI) 
#
curl "http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_prev_name&col=gd_aliases&col=gd_name_aliases&col=gd_pub_chrom_map&col=gd_pub_eg_id&col=md_eg_id&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&hgnc_dbtag=on&submit=submit" \
> HGNC_custom_download_${CURRENT_DATE}.txt

#
# Update CGD with Perl script
#
curl http://www.bbmriwiki.nl/svn/ngs_scripts/trunk/updateCGD.pl > updateCGD.pl
./updateCGD.pl \
-i CGD_${CURRENT_DATE}.txt \
-f cgd \
-sc genelist2HGNC_${CURRENT_DATE}.txt \
-cd HGNC_custom_download_${CURRENT_DATE}.txt \
-o CGD_${CURRENT_DATE}_updated.txt

#
# Tail from log from 2014-03-12:
#
2014/03/19 10:09:48 L:568 INFO> Parsed CGD file.
2014/03/19 10:09:48 L:569 INFO> ==============
2014/03/19 10:09:48 L:570 INFO> Problematic gene symbols:       0.
2014/03/19 10:09:48 L:571 INFO> ==============
2014/03/19 10:09:48 L:572 INFO> Resolved gene symbols:
2014/03/19 10:09:48 L:573 INFO>   Up-to-date HGNC gene symbols:................... 2757.
2014/03/19 10:09:48 L:574 INFO>   Gene symbols with fixed case sensitivity issue:  10.
2014/03/19 10:09:48 L:575 INFO>   Aliases / synomyms updated to HGNC gene symbols: 11.
2014/03/19 10:09:48 L:576 INFO>   ==============
2014/03/19 10:09:48 L:577 INFO>   HGNC IDs resolved via EntrezGene IDs curated by HGNC: 2716.
2014/03/19 10:09:48 L:578 INFO>   HGNC IDs resolved via EntrezGene IDs curated by NCBI: 62.
2014/03/19 10:09:48 L:579 INFO> ==============
2014/03/19 10:09:48 L:114 INFO> Finished!

#
# Display the differences with gene names.
#
awk -F '\t' 'BEGIN {OFS=FS} {if ($1 != $3) print $1,$2,$3,$5}' CGD_${CURRENT_DATE}_updated.txt
#
# Display the differences without gene names.
#
awk -F '\t' 'BEGIN {OFS=FS} {if ($1 != $3) print $1,$2,$3}' CGD_${CURRENT_DATE}_updated.txt

#HGNC_Symbol	HGNC_ID	GENE
C10orf11	HGNC:23405	C10ORF11
C10orf2	HGNC:1160	C10ORF2
C12orf57	HGNC:29521	C12ORF57
C12orf65	HGNC:26784	C12ORF65
C15orf41	HGNC:26929	C15ORF41
C19orf12	HGNC:25443	C19ORF12
C2orf71	HGNC:34383	C2ORF71
C5orf42	HGNC:25801	C5ORF42
SUGCT	HGNC:16001	C7ORF10
C8orf37	HGNC:27232	C8ORF37
C9orf72	HGNC:28337	C9ORF72
MT-CO3	HGNC:7422	COX3
MT-CYB	HGNC:7427	CYTB
MT-TH	HGNC:7487	TRNH
MT-TI	HGNC:7488	TRNI
MT-TK	HGNC:7489	TRNK
MT-TL1	HGNC:7490	TRNL1
MT-TP	HGNC:7494	TRNP
MT-TQ	HGNC:7495	TRNQ
MT-TS2	HGNC:7498	TRNS2
MT-TT	HGNC:7499	TRNT
