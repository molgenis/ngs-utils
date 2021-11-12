INPUT=/groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.merged.g.vcf.gz
OUTPUT=/groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-22.variant.calls.bed

#bedfile=/groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/targets.bed
#bedfile=/groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/ONC_3183841_+en-20_target_v1.bed
#bedfile=/groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/Exoom_target_v1plus50.chr22.big.bed
bedfile=/groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/Exoom_target_v1plus50.part.bed

# merge gvcfs
ml GATK/4.1.2.0-Java-1.8.0_144-unlimited_JCE
ml HTSlib/1.9-foss-2015b

#for b in "1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "Xp" "Xnp" "Y" "MT" "NC_001422.1"
#do
#        if [ -f "//groups/umcg-atd//tmp04//tmp//QXTR_276-Exoom_v1/run01//QXTR_276-Exoom_v1.batch-${b}.variant.calls.snpeff.vcf" ]
#        then
#                array_contains INPUTS "--variant //groups/umcg-atd//tmp04//tmp//QXTR_276-Exoom_v1/run01//QXTR_276-Exoom_v1.batch-${b}.variant.calls.snpeff.vcf" || INPUTS+=("--variant //groups/umcg-atd//tmp04//tmp//QXTR_276-Exoom_v1/run01//QXTR_276-Exoom_v1.batch-${b}.variant.calls.snpeff.vcf")
#        fi
#done

#java -jar $EBROOTGATK/gatk-package-4.1.2.0-local.jar GatherVcfs --help
#exit 0

java -jar $EBROOTGATK/gatk-package-4.1.2.0-local.jar GatherVcfs \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-2.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-3.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-4.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-5.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-6.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-7.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-8.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-9.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-10.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-11.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-12.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-13.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-14.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-15.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-16.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-17.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-18.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-19.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-20.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-21.variant.calls.g.vcf.gz \
--INPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.batch-22.variant.calls.g.vcf.gz \
--OUTPUT /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.merged.g.vcf.gz

tabix -p vcf /groups/umcg-atd/tmp04/projects/QXTR_276-Exoom_v1/run01/results/variants/gVCF/20182579_2021287_DNA110503_509403_QXTR276_S07604514.merged.g.vcf.gz


##



module load Python/3.6.3-foss-2015b
#export PYTHONPATH='/home/umcg-gvdvries/tools/:$PYTHONPATH'
#export PYTHONPATH='/home/umcg-gvdvries/test2/:$PYTHONPATH'
export PYTHONPATH='/groups/umcg-atd/tmp04/umcg-rkanninga/python_packages_gerben/tools/:$PYTHONPATH'
echo $PYTHONPATH

python gvcf2bed2.py \
-I ${INPUT} \
-O ${OUTPUT} \
-b ${bedfile} 

python ${HOME}/github/ngs-utils/bin/gvcf2bed.py \
-I ${INPUT} \
-O ${OUTPUT} \
-b ${bedfile} 


exit 0

optional arguments:
  -h, --help            show this help message and exit
  -I INPUT, --input INPUT
                        Input gVCF
  -O OUTPUT, --output OUTPUT
                        Output bed file
  -s SAMPLE, --sample SAMPLE
                        Sample name in VCF file to use. Will default to first
                        sample (alphabetically) if not supplied
  -q QUALITY, --quality QUALITY
                        Minimum genotype quality (default 20)
  -nq NON_VARIANT_QUALITY, --non-variant-quality NON_VARIANT_QUALITY
                        Minimum genotype quality for non-variant records
                        (default 20)
  -b, --bedgraph        Output in bedgraph mode