#!/bin/bash
#SBATCH --job-name=GavinStandAlone_VULIN
#SBATCH --output=GavinStandAlone_VULIN.out
#SBATCH --error=GavinStandAlone_VULIN.err
#SBATCH --time=12:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 11gb
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=30L

set -e
set -u

#
##Step 1, prepare
#
tmpName="tmp04"
host=$(hostname)
echo "HOST is $host"

if [[ "${host}" == "zinc-finger.gcc.rug.nl" || "${host}" == *"zf-compute"* ]]
then
	echo "tmpie05"
	tmpName="tmp05"
elif [[ "${host}" == "leucine-zipper.gcc.rug.nl" || "${host}" == *"lc-compute"* ]]
then
	echo "tmpie06"

	tmpName="tmp06"
fi
SKIPINDEL="FALSE"
SKIPSNP="FALSE"

GAVINHOME=/groups/umcg-gd/${tmpName}/GavinStandAlone/
echo "GAVINHOME: $GAVINHOME"
vcfFile=VULHIERJEVCFFILEIN
bas=$(basename ${vcfFile})
vcfName=${bas%%.*}

tmpFolder=${GAVINHOME}/tmp/${vcfName}/
jobsFolder=${GAVINHOME}/jobs/${vcfName}/

snpeff=${tmpFolder}/${vcfName}.snpeff.vcf
exac=${tmpFolder}/${vcfName}.ExAC_Annotated.vcf
gonl=${tmpFolder}/${vcfName}.ExAC_GoNL_Annotated.vcf
cadd=${tmpFolder}/${vcfName}.ExAC_GoNL_CADD_Annotated.vcf
snps=${tmpFolder}/${vcfName}_snps.vcf
indels=${tmpFolder}/${vcfName}_indels.vcf
snpsFiltered=${snps%.*}_filtered.vcf
indelsFiltered=${indels%.*}_filtered.vcf

sleep 2

mkdir -p ${tmpFolder}
mkdir -p ${jobsFolder}

awk '{if ($1 == "#CHROM"){print $0}}' ${vcfFile} | awk '{for(i=10;i<=NF;++i)print $i}' > ${tmpFolder}/externalSampleIDs.txt

size=$(cat ${tmpFolder}/externalSampleIDs.txt | wc -l)

if [ ${size} == 1 ]
then
	echo "${vcfName}" > ${tmpFolder}/externalSampleIDs_fix.txt 
else
	cat ${tmpFolder}/externalSampleIDs.txt > ${tmpFolder}/externalSampleIDs_fix.txt
fi
#
##Step 2, annotating with SnpEff
#
mv $vcfFile ${GAVINHOME}/input/processing/
echo "moving $vcfFile to ${GAVINHOME}/input/processing/"

if [ ! -f ${tmpFolder}/snpeff.sh.finished ]
then

	module load snpEff/4.3-Java-1.7.0_80
	module list

    	#Run snpEff
        java -XX:ParallelGCThreads=4 -Xmx4g -jar \
        $EBROOTSNPEFF/snpEff.jar \
        -v hg19 \
        -noStats \
        -noLog \
        -lof \
	-canon \
        -ud 0 \
        -c $EBROOTSNPEFF/snpEff.config \
        ${GAVINHOME}/input/processing/${bas} \
        > ${snpeff}~

        mv ${snpeff}~ ${snpeff}
        touch ${tmpFolder}/snpeff.sh.finished
else
    	echo "snpeff skipped"
fi

if [ ! -f ${tmpFolder}/cmdlineannotator.sh.finished ]
then

	ml CmdLineAnnotator

    	echo "exac annotation"
                java -Xmx10g -jar ${EBROOTCMDLINEANNOTATOR}/CmdLineAnnotator-1.21.1.jar \
                -a exac \
                -s //apps//data//ExAC/release0.3/ExAC.r0.3.sites.vep.vcf.gz \
                -i ${snpeff} \
                -o ${exac}~

                mv ${exac}~ ${exac}
                echo "mv ${exac}~ ${exac}"

        echo "exac done"
        echo "gonl annotation"

                java -Xmx10g -jar ${EBROOTCMDLINEANNOTATOR}/CmdLineAnnotator-1.21.1.jar \
                -a gonl \
                -s //apps//data//gonl/release5_noContam_noChildren_with_AN_AC_GTC_stripped/ \
                -i ${exac} \
                -o ${gonl}~

                mv ${gonl}~ ${gonl}

                echo "mv ${gonl}~ ${gonl}"

        echo "gonl done"
        echo "cadd annotation"

                java -Xmx10g -jar ${EBROOTCMDLINEANNOTATOR}/CmdLineAnnotator-1.21.1.jar \
                -a cadd \
                -s //apps//data//CADD/whole_genome_SNVs.tsv.gz \
                -i ${gonl} \
                -o ${cadd}~

                mv ${cadd}~ ${cadd}
                echo "mv ${cadd}~ ${cadd}"
        echo "cadd done"
        perl -i -wpe 'my @t = split("\t",$_);$t[7] =~ s/ /_/g if ($t[7]);$t[7] =~ s/;$//g if ($t[7]);$_ = join("\t",@t)' ${cadd}

        touch ${tmpFolder}/cmdlineannotator.sh.finished
else
    	echo "cmdlineannotator skipped"
fi


if [ ! -f ${tmpFolder}/splitsnp.sh.finished ]
then
	
	module load GATK/3.6-Java-1.8.0_74

    	java -XX:ParallelGCThreads=2 -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
        -R //apps//data//1000G/phase1/human_g1k_v37_phiX.fasta \
        -T SelectVariants \
        --variant ${cadd} \
        -o ${snps} \
        -L ${cadd} \
        --selectTypeToExclude INDEL

    	perl -pi -e 's|ID=GC,Number=1,Type=Integer|ID=GC,Number=1,Type=Float|g' ${snps}
        perl -pi -e 's|MQ0Fraction,Number=1,Type=Integer|MQ0Fraction,Number=1,Type=Float|g' ${snps}
        perl -pi -e 's|MQ0,Number=1,Type=Integer|MQ0,Number=1,Type=Float|g' ${snps}
	
	awk '{if ($1 !~ /^#/){print $0}}' ${snps} > ${snps}.count	
	size=$(cat ${snps}.count | wc -l)
	if [ $size == 0 ]
	then
		SKIPSNP="TRUE"
		touch ${tmpFolder}/skipsnp.true
	fi	
        touch ${tmpFolder}/splitsnp.sh.finished
else
    	echo "splitsnp already finished"
fi

if [ ! -f ${tmpFolder}/splitindel.sh.finished ]
then
	module load GATK/3.6-Java-1.8.0_74

    	java -XX:ParallelGCThreads=2 -Xmx4g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
        -R //apps//data//1000G/phase1/human_g1k_v37_phiX.fasta \
        -T SelectVariants \
        --variant ${cadd} \
        -o ${indels} \
        -L ${cadd} \
        --selectTypeToInclude INDEL

	perl -pi -e 's|ID=GC,Number=1,Type=Integer|ID=GC,Number=1,Type=Float|g' ${indels}     
        perl -pi -e 's|MQ0Fraction,Number=1,Type=Integer|MQ0Fraction,Number=1,Type=Float|g' ${indels}       
	perl -pi -e 's|INFO=\<ID=MQ0,Number=1,Type=Integer|INFO=\<ID=MQ0,Number=1,Type=Float|g' ${indels}

	awk '{if ($1 !~ /^#/){print $0}}' ${indels} > ${indels}.count
	size=$(cat ${indels}.count | wc -l)

        if [ $size == 0 ]
        then
                SKIPINDEL="TRUE"
		touch ${tmpFolder}/skipindel.true
        fi

        touch ${tmpFolder}/splitindel.sh.finished
else
    	echo "splitindel finished"
fi

if [ -f ${tmpFolder}/skipindel.true ] 
then 
	SKIPINDEL="TRUE"
fi
if [ -f ${tmpFolder}/skipsnp.true ] 
then 
	SKIPSNP="TRUE"
fi


if [ ! -f ${tmpFolder}/indelsnpfiltration.sh.finished ]
then

	ml GATK/3.6-Java-1.8.0_74

	inputSNPwithoutextension=${snps%.*}
	inputIndelwithoutextension=${indels%.*}

	if [ "${SKIPINDEL}" == "FALSE" ]
	then

    		java -XX:ParallelGCThreads=4 -Xmx8g -Xms6g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
	        -T VariantFiltration \
	        -R //apps//data//1000G/phase1/human_g1k_v37_phiX.fasta \
	        -o ${indelsFiltered} \
	        --variant ${indels} \
	        --filterExpression "QD < 2.0" \
	        --filterName "filterQD" \
	        --filterExpression "FS > 200.0" \
	        --filterName "filterFS" \
	        --filterExpression "ReadPosRankSum < -20.0" \
	        --filterName "filterReadPosRankSum"
	
	fi
	
	if [ "${SKIPSNP}" == "FALSE" ]
        then
	        java -XX:ParallelGCThreads=4 -Xmx8g -Xms6g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
	        -T VariantFiltration \
	        -R //apps//data//1000G/phase1/human_g1k_v37_phiX.fasta \
	        -o ${snpsFiltered} \
	        --variant ${snps}   \
	        --filterExpression "QD < 2.0" \
	        --filterName "filterQD" \
	        --filterExpression "MQ < 25.0" \
	        --filterName "filterMQ" \
		--filterExpression "FS > 60.0" \
	        --filterName "filterFS" \
	        --filterExpression "MQRankSum < -12.5" \
	        --filterName "filterMQRankSum" \
	        --filterExpression "ReadPosRankSum < -8.0" \
	        --filterName "filterReadPosRankSum"
	
	fi

        touch ${tmpFolder}/indelsnpfiltration.sh.finished

else
    	echo "indelsnpfiltration skipped"
fi

if [ ! -f ${tmpFolder}/mergeindelsnp.sh.finished ]
then

	if [[ "${SKIPSNP}" == "FALSE" && "${SKIPINDEL}" == "FALSE"  ]]
	then

		module load GATK/3.6-Java-1.8.0_74

		INPUTS=()

		while read line
		do
		  	INPUTS+=("$line")
		done<${tmpFolder}/externalSampleIDs.txt

	    	for externalID in "${INPUTS[@]}"
	        do
	          	java -Xmx2g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
	                -R //apps//data//1000G/phase1/human_g1k_v37_phiX.fasta \
	                -T CombineVariants \
	                --variant $snpsFiltered \
	                --variant $indelsFiltered \
	                --genotypemergeoption UNSORTED \
	                -o ${tmpFolder}/${externalID}.final.vcf

	        done
	elif [ "${SKIPINDEL}" == "TRUE" ]
        then
            	INPUTS=()

                while read line
                do
                  	INPUTS+=("$line")
                done<${tmpFolder}/externalSampleIDs.txt

                for externalID in "${INPUTS[@]}"
                do
                  	cat ${snpsFiltered} > ${tmpFolder}/${externalID}.final.vcf
                done
        elif [ "${SKIPSNP}" == "TRUE" ]
        then
            	INPUTS=()

                while read line
                do
                        INPUTS+=("$line")
                done<${tmpFolder}/externalSampleIDs.txt

            	for externalID in "${INPUTS[@]}"
                do
                  	cat $indelsFiltered > ${tmpFolder}/${externalID}.final.vcf
                done
        fi
	touch ${tmpFolder}/mergeindelsnp.sh.finished
else
    	echo "mergeindelsnp skipped"
fi


if [ ! -f ${tmpFolder}/selectvariants.finished ]
then
	module load GATK/3.6-Java-1.8.0_74

        INPUTS=()

        while read line
        do
                INPUTS+=("$line")
        done<${tmpFolder}/externalSampleIDs.txt
    	
	for externalID in "${INPUTS[@]}"
        do
          	java -Xmx10g -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
                -R //apps//data//1000G/phase1/human_g1k_v37_phiX.fasta \
                -V ${tmpFolder}/${externalID}.final.vcf \
                -T SelectVariants \
                -o ${tmpFolder}/${externalID}.splitted.final.vcf \
                -sn $externalID
	
        done
 	if [[ ${size} -ne 1 ]]
	then
	    	echo "skippy"
	else
    		cat ${tmpFolder}/externalSampleIDs.txt > ${tmpFolder}/externalSampleIDs_fix.txt

		for i in $(ls ${tmpFolder}/*.splitted.final.vcf)
		do
			baseNameSplitted=$(basename $i)
			vcfNameSplitted=${baseNameSplitted%%.*}
			externalUpdatedID=$(ls ${tmpFolder}/*${vcfNameSplitted}*.snpeff.vcf)
	
			basenameUpdated=$(basename $externalUpdatedID)
			vcfNameUpdated=${basenameUpdated%%.*}
			echo "ombatterijen van splitted.final.vcf"	
		
			if [ ${i} == ${tmpFolder}/${vcfNameUpdated}.splitted.final.vcf ]
			then
				echo "files are the same, no moving of splitted.vcf"
			else
				mv ${i} ${tmpFolder}/${vcfNameUpdated}.splitted.final.vcf
				mv ${tmpFolder}/${vcfNameSplitted}.final.vcf ${tmpFolder}/${vcfNameUpdated}.final.vcf
				mv ${i}.idx ${tmpFolder}/${vcfNameUpdated}.splitted.final.vcf.idx
				mv ${tmpFolder}/${vcfNameSplitted}.final.vcf.idx ${tmpFolder}/${vcfNameUpdated}.final.vcf.idx
				echo "moving ${tmpFolder}/${vcfNameSplitted}.final.vcf ${tmpFolder}/${vcfNameUpdated}.final.vcf DONE"

			fi
		
		done
	fi
	touch ${tmpFolder}/selectvariants.finished
else
    	echo "selectvariants skipped"
fi
INPUTS=()

while read line
do
	INPUTS+=("$line")
done<${tmpFolder}/externalSampleIDs_fix.txt

for externalID in "${INPUTS[@]}"
do
echo -e "#!/bin/bash
#SBATCH --job-name=${externalID}_Gavin
#SBATCH --output=${externalID}_Gavin.out
#SBATCH --error=${externalID}_Gavin.err
#SBATCH --time=65:59:00
#SBATCH --cpus-per-task 1
#SBATCH --mem 11gb
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
	
set -e
set -u
	module load HTSlib/1.3.2-foss-2015b

	module load Gavin-ToolPack/1.0-Java-1.8.0_74
	input=${tmpFolder}/${externalID}.splitted.final.vcf
	sampleID=${externalID}
	perl -pi -e 's|ID=GC,Number=1,Type=Integer|ID=GC,Number=1,Type=Float|g' \${input}
        perl -pi -e 's|MQ0Fraction,Number=1,Type=Integer|MQ0Fraction,Number=1,Type=Float|g' \${input}
        perl -pi -e 's|INFO=\<ID=MQ0,Number=1,Type=Integer|INFO=\<ID=MQ0,Number=1,Type=Float|g' \${input}

    	module list

	java -Xmx10g -jar \${EBROOTGAVINMINTOOLPACK}/GAVIN-APP-1.0.jar \\
        -i \${input} \\
        -o \${input}.firstpass~ \\
        -m CREATEFILEFORCADD \\
        -r \\
	-a \${input}.toCadd.tsv.tmp \\
        -c //apps//data//GAVIN//clinvar.patho.fix.11oct2016.vcf.gz \\
        -d //apps//data//GAVIN//CGD_11oct2016.txt.gz \\
        -f //apps//data//GAVIN//FDR_allGenes_r1.0.tsv \\
        -g //apps//data//GAVIN//GAVIN_calibrations_r0.3.tsv

        echo \${input}.toCadd.tsv.tmp
        awk '{if (\$1 != \"NC_001422.1\"){print \$0}}' \${input}.toCadd.tsv.tmp >  \${input}.toCadd.tsv
        mv \${input}.firstpass~ \${input}.firstpass

        echo \"GAVIN round 1 is finished, uploading to CADD...\"

        echo \"starting to get CADD annotations locally for \${input}.toCadd.tsv\"
        if [ ! -f \${input}.fromCadd.tsv.gz ]
        then
            	bgzip -c \${input}.toCadd.tsv > \${input}.toCadd.tsv.gz
                tabix -p vcf \${input}.toCadd.tsv.gz

                score.sh \${input}.toCadd.tsv.gz \${input}.fromCadd.tsv.gz
		gzip -d \${input}.fromCadd.tsv.gz
        else
            	echo \"CADD already done, skipped\"
        fi

	java -Xmx10g -jar \${EBROOTGAVINMINTOOLPACK}/GAVIN-APP-1.0.jar \\
        -i \${input} \\
        -o \${input}.final.vcf \\
        -m ANALYSIS \\
        -r \\
	-a \${input}.fromCadd.tsv \\
        -c //apps//data//GAVIN//clinvar.patho.fix.11oct2016.vcf.gz \\
        -d //apps//data//GAVIN//CGD_11oct2016.txt.gz \\
        -f //apps//data//GAVIN//FDR_allGenes_r1.0.tsv \\
        -g //apps//data//GAVIN//GAVIN_calibrations_r0.3.tsv

        echo \"Merging \${input} and \${input}.final.vcf\"	

	java -jar -Xmx10g \${EBROOTGAVINMINTOOLPACK}/MergeBackTool-0.2.jar \\
        -i \${input} \\
        -r \\
	-v \${input}.final.vcf \\
        -o \${input}.GAVIN.RVCF.final.mergedWithOriginal.vcf


        #gavinSplitRlvToolJar
        java -jar -Xmx10g \${EBROOTGAVINMINTOOLPACK}/SplitRlvTool-0.2.jar \\
        -i \${input}.GAVIN.RVCF.final.mergedWithOriginal.vcf \\
        -r \\
	-o \${input}.GAVIN.RVCF.final.mergedWithOriginal.rlv.vcf

        perl -pi -e 's|INFO=<ID=EXAC_AF,Number=.,Type=String|INFO=<ID=EXAC_AF,Number=.,Type=Float|' \${input}.GAVIN.RVCF.final.mergedWithOriginal.rlv.vcf
        perl -pi -e 's|INFO=<ID=EXAC_AC_HOM,Number=.,Type=String|INFO=<ID=EXAC_AC_HOM,Number=.,Type=Integer|' \${input}.GAVIN.RVCF.final.mergedWithOriginal.rlv.vcf
        perl -pi -e 's|INFO=<ID=EXAC_AC_HET,Number=.,Type=String|INFO=<ID=EXAC_AC_HET,Number=.,Type=Integer|' \${input}.GAVIN.RVCF.final.mergedWithOriginal.rlv.vcf
	mv \${input}.GAVIN.RVCF.final.mergedWithOriginal.rlv.vcf ../output/${externalID}.GAVIN.rlv.vcf
	mv ../input/processing/${bas} ../input/done
	mv ../tmpGavin/${externalID}_Gavin.{err,out} ../tmp/${externalID}/
        
	fi

        touch ${tmpFolder}/gavin.${externalID}.sh.finished

"> ${jobsFolder}/gavin.${externalID}.sh

done

for externalID in "${INPUTS[@]}"
do
  	if [ ! -f tmpGavin/gavin.${externalID}.sh.finished ]
        then
		sbatch  ${jobsFolder}/gavin.${externalID}.sh ${tmpFolder}/${externalID}.splitted.final.vcf $externalID

	else
    		echo "gavin.${externalID} already finished, skipped"
	fi

done
