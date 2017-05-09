region=$1
bam=$2
ml GATK/3.7-Java-1.8.0_74 
gender=""
if [[ $region == *"X"* || $region == *"x"* ]]
then
	if [ -z "$3" ]
  	then
    		echo "No gender supplied, please fill in a gender (Male or Female) when using chromosome X region"
	else

		gender=$3
	fi
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R human_b37_20.fasta -I recalibrated.bam -o hc_variants.vcf -L 20:10255630-10255840 -bamout bamout.bam

java -XX:ParallelGCThreads=2 -Xmx12g -jar \
${EBROOTGATK}/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R //apps//data//1000G/phase1/human_g1k_v37_phiX.fasta \
-I ${bam} \
-dontUseSoftClippedBases \
--dbsnp //apps//data//dbSNP//dbsnp_137.b37.vcf \
-bamout ${bam}.${region}.bam
-o "${bam}".vcf \
-L "${region}"
--emitRefConfidence GVCF \
-ploidy 2
