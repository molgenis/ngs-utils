#!/bin/bash
#SBATCH --job-name=ConcordanceCheck
#SBATCH --output=ConcordanceCheck.out
#SBATCH --error=ConcordanceCheck.err
#SBATCH --partition=devel
#SBATCH --time=43:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 6gb
#SBATCH --nodes 1
#SBATCH --open-mode=append

set -e
set -u

### Parameters
. ConcordanceCheckConfig

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
	--recode-vcf \
	--ped ${s}.concordance.ped \
	--map ${arrayMapFile} \
	--out ${s}.concordance

	##Rename plink.vcf to sample.vcf
	mv ${s}.concordance.vcf ${s}.genotypeArray.vcf

	##Replace chr23 and 24 with X and Y
    	perl -pi -e 's/^23/X/' ${s}.genotypeArray.vcf
 	perl -pi -e 's/^24/Y/' ${s}.genotypeArray.vcf
	
	##Create binary ped (.bed) and make tab-delimited .fasta file for all genotypes
	sed -e 's/chr//' ${sample}.genotypeArray.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2}}' \
	> ${s}.genotypeArray.bed

	###########################################################################################

	### Compare Array data with NGS vcf-File using GATK GenotypeConcordance


	#Removing small indels from NGS VCF
	
	module load GATK/3.6-Java-1.8.0_74
	 java -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
        -T SelectVariants \
   	-R ${r} \
   	-V ${v} \
   	-o ${n}.onlySNPs.vcf \
   	-selectType SNP \
	
	
	### Comparing VCF From NGS with Array VCF
	module load GATK/3.6-Java-1.8.0_74 

	java -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
	-T GenotypeConcordance \
	-R ${r} \
	-eval  ${n}.onlySNPs.vcf  \
	-comp ${s}.genotypeArray.vcf \
	-o ${s}.GATK.VCF.Concordance.output.grp 


	### Compare Array data with raw NGS data using verifyBamID
	
	module load verifyBamID/1.1.2-goolf-1.7.20
	
	verifyBamID \
	--bam ${b} \
	--vcf ${s}.genotypeArray.vcf \
	--out ${s}.VerifyBamID \
	--verbose \
	--ignoreRG \
	--chip-mix
	
	exit 1


fi
