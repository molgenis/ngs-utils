import re
entities = ['Adapter_content','AddOrReplaceReadGroups','Alignment_summary_metrics','Alts','AnalyseCovariates','Arguments','BaseQualityScoreRecalibration','BQSR','BQSR_after_grp','BQSR_before_grp','CMMetrics','CombineBedFiles','Comp_overlap','Contigs','Count_variants','CRMetrics','DepthRG','DepthSM','Duplication_metrics','ENA','FastQC','Filters','Flagstat','Formats','GATKSplitNTrim','GenotypeGvcf','GenotypeHarmonizer','Gvcf_block','Hisat','Indel_length_histogram','Indel_summary','IndelRealignmentKnown','Individuals','Info','InHouse','Insert_size_metrics_class','Insert_size_metrics_histogram','Kmer_content','MarkDuplicates','MarkDuplicates_histogram','Md5sums','MergeBamFiles','MergeGvcf','Multiallelic_summary','Overrepresented_sequences','Per_base_N_content','Per_base_sequence_content','Per_base_sequence_quality','Per_sequence_GC_content','Per_sequence_quality_scores','Per_tile_sequence_quality','Quality_by_cycle_metrics','Quality_distribution_metrics','Quantized','Recal_table_0','Recal_table_2','Recal_table_2','Rnaseq_metrics_class','Rnaseq_metrics_histogram','Samples','SamToFilteredBam','SelfRG','SelfSM','Sequence_duplication_levels','Sequence_length_distribution','Skipped_marker','Snp_info','Ti_tv_variant_evaluator','Tissue','Validation_report','Variant_class_count','Variant_summary','VariantCaller','VariantEval','Vcf','VerifyBamID','VerifyBamID_individual','SortBam']
total_entities = len(entities)
# entities in variables, won't be found by this code
entities_done = ['DepthSM','SelfSM','Variant_summary','Validation_report','Ti_vv_variant_evaluator',
                 'Multiallelic_summary','IndelLength_histogram','Count_variants','Comp_overlap','Indel_summary',
                 'Variant_class_count','BQSR_before_grp','BQSR_after_grp','DepthRG','DepthSM','Ti_tv_variant_evaluator',
                 'CMMetrics','BaseQualityScoreRecalibration','GenotypeGvcf','SelfRG','Indel_length_histogram']
entities = list(set(entities)-set(entities_done))
database_filler = open('RNAseqParser/parse_output.py')

packages = []
for line in database_filler:
    line = line.strip()
    if 'connection.add_entity_row' in line and not 'package+entity' in line:
        try: 
            print(line)
            package = re.search('connection.add_entity_row\(package\+\'(\w+)',line).group(1)
            packages.append(package)
        except IndexError:
            packages.append(line)
print((str(len(set(packages))+len(entities_done))+'/'+str(total_entities)))
print((set(entities) - set(packages)))