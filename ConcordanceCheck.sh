#!/bin/bash
#SBATCH --job-name=ConcordanceCheck
#SBATCH --output=ConcordanceCheck.out
#SBATCH --error=ConcordanceCheck.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 15gb
#SBATCH --nodes 1
#SBATCH --open-mode=append

set -e
set -u

### Parameters
. $PWD/ConcordanceCheckConfig.cfg

###################################################################################

###Start protocol


if test ! -e ${f};
then
	echo "name, step, nSNPs, PercDbSNP, Ti/Tv_known, Ti/Tv_Novel, All_comp_het_called_het, Known_comp_het_called_het, Non-Ref_Sensitivity, Non-Ref_discrepancy, Overall_concordance" > ${sampleConcordanceFile}
	echo "[1] NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA" >> ${sampleConcordanceFile} 
else
	#Check finalReport on "missing" alleles. Also, see if we can fix missing alleles somewhere in GenomeStudio
	awk '{ if ($3 != "-" || $4 != "-") print $0};' ${f} \
	> ${s}_FinalReport.txt.tmp
	
	#Check finalreport on "D"alleles.
	awk '{ if ($3 != "D" || $4 != "D") print $0};' ${s}_FinalReport.txt.tmp \
	> ${s}_FinalReport_2.txt.tmp

	#Push sample belonging to family "1" into list.txt
	echo 1 ${s} > ${familyList}
	
	#########################################################################
	#########################################################################
	
	module load ngs-utils/16.06.1 

	##Create .fam, .lgen and .map file from sample_report.txt
	sed -e '1,10d' ${s}_FinalReport_2.txt.tmp | awk '{print "1",$2,"0","0","0","1"}' | uniq > ${s}.concordance.fam
	sed -e '1,10d' ${s}_FinalReport_2.txt.tmp | awk '{print "1",$2,$1,$3,$4}' | awk -f ${toolDir}RecodeFRToZero.awk > ${s}.concordance.lgen
	sed -e '1,10d' ${s}_FinalReport_2.txt.tmp | awk '{print $6,$1,"0",$7}' OFS="\t" | sort -k1n -k4n | uniq > ${arrayTmpMap}
	grep -P '^[123456789]' ${arrayTmpMap} | sort -k1n -k4n > ${arrayMapFile}
	grep -P '^[X]\s' ${arrayTmpMap} | sort -k4n >> ${arrayMapFile}
	grep -P '^[Y]\s' ${arrayTmpMap} | sort -k4n >> ${arrayMapFile}

	#####################################
	##Create .bed and other files (keep sample from sample_list.txt).
	
	##Create .bed and other files (keep sample from sample_list.txt).
	
	module load PLINK/1.07-x86_64
	module list
	
	plink \
	--lfile ${s}.concordance \
	--recode \
	--noweb \
	--out ${s}.concordance \
	--keep ${familyList}

	module unload plink
	module load plink/1.9-foss-2015b 
	module list
	
	##Create genotype VCF for sample
	plink \
	--recode vcf-iid \
	--ped ${s}.concordance.ped \
	--map ${arrayMapFile} \
	--out ${s}.concordance

	##Rename plink.vcf to sample.vcf
	mv ${s}.concordance.vcf ${s}.genotypeArray.vcf

	##Replace chr23 and 24 with X and Y
  	perl -pi -e 's/^23/X/' ${s}.genotypeArray.vcf
	perl -pi -e 's/^24/Y/' ${s}.genotypeArray.vcf
	
	##Create binary ped (.bed) and make tab-delimited .fasta file for all genotypes
	sed -e 's/chr//' ${s}.genotypeArray.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2}}' \
	> ${s}.genotypeArray.bed

	
	#Remove SNP`s from array which are not in a exon with the exon bedfile
	
	module load BEDTools/2.25.0-foss-2015b
	bedtools intersect -a ${s}.genotypeArray.vcf -b /apps/data/UMCG/Ensembl.GRCh37.75-AllExons_+20_-20bp_AgilentV5_AllExon_SSID_CGD_2015-11-26/Ensembl.GRCh37.75-AllExons_+20_-20bp_AgilentV5_AllExon_SSID_CGD_2015-11-26.bed -header  >${s}.genotypeArray.ExonFiltered.vcf


	#Remove SNP's from array which are called homozygous reference
	awk '{ if ($10!= "0/0") print $0};' ${s}.genotypeArray.ExonFiltered.vcf \
	> ${s}.genotypeArray.ExonFiltered.HomozygousRefRemoved.vcf 

	sleep 3m
		
	#Count how much SNP's are in original VCF and how much proceed for Concordance
	wc -l ${s}.genotypeArray.vcf > originalSNPs.txt
	wc -l ${s}.genotypeArray.ExonFiltered.HomozygousRefRemoved.vcf > SNPswichproceedtoConcordance.txt

        #Change Array VCF to same name as NGS VCF
        awk '{OFS="\t"}{if ($0 ~ "#CHROM" ){ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'$n'"} else {print $0}}' ${s}.genotypeArray.ExonFiltered.HomozygousRefRemoved.vcf  > ${n}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf


        #Making Array VCF index

        module load tabix/0.2.6-foss-2015b
	bgzip -c ${n}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf > ${n}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf.gz
        tabix -p vcf ${n}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf.gz

	#Uncompress Files

	#bgzip -d ${n}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf.gz.tbi > ${n}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf.tbi

	###########################################################################################
	#Updating NGS VCF File

        #Remove reference calls from NGS VCF

 #       module load VCFtools/0.1.12b-foss-2015b-Perl-5.20.2-bare
 #       vcf-subset -e ${v} > ${n}.removed_reference_calls.vcf

	#Removing small indels from NGS VCF
	
	module load GATK/3.6-Java-1.8.0_74
	java -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
        -T SelectVariants \
  	-R ${r} \
   	-V ${v} \
  	-o ${n}.onlySNPs.FINAL.vcf \
   	-selectType SNP

	
	########################################################################################################

	### Compare Array data with NGS vcf-File using GATK GenotypeConcordance
	
	### Comparing VCF From NGS with Array VCF
	module load GATK/3.6-Java-1.8.0_74 

	java -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
	-T GenotypeConcordance \
	-R ${r} \
	-eval ${n}.onlySNPs.FINAL.vcf \
	-comp ${n}.genotypeArray.ExonFiltered.HomozygousRefRemoved.FINAL.vcf \
	-o ${n}.GATK.VCF.Concordance.output.grp 


	### Compare Array data with raw NGS data using verifyBamID
	
#	module load verifyBamID/1.1.2-goolf-1.7.20
#	
#	verifyBamID \
#	--bam ${b} \
#	--vcf ${s}.genotypeArray.ExonFiltered.vcf \
#	--out ${s}.VerifyBamID \
#	--smID ${s} \
#	--verbose \
#	--ignoreRG \
#	--best \
#	--chip-mix

##Put Relevant outputdata to output folder

	mkdir Output_${n}_Concordance
	cp ${n}.GATK.VCF.Concordance.output.grp Output_${n}_Concordance
	cp SNPswichproceedtoConcordance.txt Output_${n}_Concordance
	cp originalSNPs.txt Output_${n}_Concordance
	 
	
	
fi
