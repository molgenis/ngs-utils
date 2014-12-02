#
# ======================================================
# How to update gene symbols
# referenced in the Clinical Genomics Database (CGD)
# to official gene symbols as curated by the HGNC
#
# Last update: 2014 October 08
# ======================================================
#

#CGD_DATE=$(date "+%Y-%m-%d")
CGD_DATE=2014-10-08

#
# Fetch new CGD from research.nhgri.nih.gov/CGD
#
curl http://research.nhgri.nih.gov/CGD/download/txt/CGD.txt.gz > CGD_${CGD_DATE}.txt.gz
gunzip CGD_${CGD_DATE}.txt.gz

#
# Create genesymbol2HGNC.txt by uploading gene symbols from CGD to HGNC gene symbol checker
#
# Go to http://www.genenames.org/cgi-bin/symbol_checker
# Paste the gene symbols from CGD_${CGD_DATE}.txt
# Save as genelist2HGNC_${CGD_DATE}.txt 
#
# Fetch "custom download" with Entrez Gene IDs
#
# Go to http://www.genenames.org/cgi-bin/download
# Select columns
#   HGNC ID
#   Approved Symbol
#   Status
#
curl "http://www.genenames.org/cgi-bin/download?col=gd_hgnc_id&col=gd_app_sym&col=gd_status&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&submit=submit" \
> HGNC_custom_download_${CGD_DATE}.txt

#
# Swap first 2 columns to make sure the first column contains the HGNC ID, 
# which is used for the join later on.
#
awk -F '\t' 'BEGIN {OFS=FS} {temp = $1; $1 = $2; $2 = temp; print;}' \
   CGD_${CGD_DATE}.txt \
 > CGD_${CGD_DATE}.txt.swapped
mv CGD_${CGD_DATE}.txt{.swapped,}

#
# Update CGD by joining on HGNC ID on the commandline.
#
join -t $'\t' -a 1 <(sort -k 1,1 CGD_${CGD_DATE}.txt) <(sort -k 1,1 HGNC_custom_download_${CGD_DATE}.txt) \
 | sort -k1,1n \
 > CGD_${CGD_DATE}_updated.txt

#
# Display the differences.
#
awk -F '\t' 'BEGIN {OFS=FS} {if ($2 != $13) print $1,$2,$13}' CGD_${CGD_DATE}_updated.txt

#HGNC ID	GENE	Approved Symbol
#1160	C10ORF2	C10orf2
#7422	COX3	MT-CO3
#7427	CYTB	MT-CYB
#7487	TRNH	MT-TH
#7488	TRNI	MT-TI
#7489	TRNK	MT-TK
#7490	TRNL1	MT-TL1
#7494	TRNP	MT-TP
#7495	TRNQ	MT-TQ
#7498	TRNS2	MT-TS2
#7499	TRNT	MT-TT
#23405	C10ORF11	C10orf11
#25443	C19ORF12	C19orf12
#25801	C5ORF42	C5orf42
#26784	C12ORF65	C12orf65
#26929	C15ORF41	C15orf41
#27232	C8ORF37	C8orf37
#28337	C9ORF72	C9orf72
#29521	C12ORF57	C12orf57
#34383	C2ORF71	C2orf71

#
# Create list of HGNC Gene Symbols present in CGD.
#
awk -F '\t' 'BEGIN {OFS=FS} {print $13}' CGD_${CGD_DATE}_updated.txt | sort > CGD_${CGD_DATE}_updated.genelist


