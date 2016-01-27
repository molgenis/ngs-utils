'''
Created on Jul 4, 2015

@author: Niek

TODO: Make sure all non-rnaseqtool tables also get duplicate value prevention
'''
from RNAseqParser import molgenis_wrapper
import warnings
import os
import re
import glob
import datetime
import zipfile
import configparser
import requests
config = configparser.RawConfigParser()
config.read(r'RNAseqParser/CONFIG')
if __name__ == "__main__":
    config.read(r'CONFIG')
def configSectionMap(section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        dict1[option] = config.get(section, option)
        if dict1[option] == -1:
            print(("skip: %s" % option))
    return dict1
analysis_id = configSectionMap("settings")['analysis_id']
analysis_description = configSectionMap("settings")['analysis_description']
def add_multiple_rows(entity, data,connection):
    '''Helper function for adding multiple rows at the same time. Takes the max_rows from the CONFIG file
       and a list of data to be added to entity, adds them all, and returns a list of all the IDs
    
    Args:
        entity (string): Entity to add it to
        data (list): list of data to add to entity
        connection (obj): Connection object
    Returns:
        list of added IDs
    '''
    def chunker(seq, size):
        return (seq[pos:pos + size] for pos in range(0, len(seq), size))
    added_ids = []
    for data in chunker(data,int(configSectionMap("settings")['max_rows'])):
        # merge the lists of returned IDs
        added_ids += connection.add(entity,data=data)
    return added_ids
def parse_ena(ena_file,connection,package,project):
    '''finished'''
    print ('Start ENA')
    columns_to_change = ['fastq_aspera','fastq_bytes','fastq_ftp','fastq_md5']
    
    with open(ena_file,'r', encoding="utf-8") as ena_meta:
        header = ena_meta.readline().strip().split('\t')
        to_add = []
        count = 0
        for line in ena_meta:
            count += 1
            column_data = line.strip().split('\t')
            library_layout = column_data[12]
            data = {}
            for column, value in zip(header, column_data):
                # columns to change are those that are different depending on if it's first pair or second pair
                if len(column_data) < 30:
                    data['fastq_files_available'] = False
                    data['correct_library_layout'] = "UNKNOWN"
                else: 
                    data['fastq_files_available'] = True
                    if column in columns_to_change:
                        if ';' in column_data[29] and '_2.fastq.gz' in column_data[29] and '_1.fastq.gz' in column_data[29]:
                            # the paired columns are not always in the same order in the ENA file
                            data['correct_library_layout'] = 'PAIRED'
                            if column_data[29].split(';')[0].endswith('_1.fastq.gz'):
                                first_of_pair = 0
                                second_of_pair = 1
                            elif column_data[29].split(';')[0].endswith('_2.fastq.gz'):
                                first_of_pair = 1
                                second_of_pair = 0
                            if ';' in value:
                                split_value = value.split(';')
                                value = split_value[first_of_pair]
                                data[column+'_2'] = split_value[second_of_pair]
                        else:
                            data['correct_library_layout'] = 'SINGLE'
    
                        column += '_1'
                data[column] = value
            to_add.append(data)
        added_ids = add_multiple_rows(entity=package+'ENA',data=to_add,connection=connection)
        for i, data in enumerate(to_add):
            connection.update_entity_rows(package+'Samples', data={'verifyBamID':added_ids[i]}, row_id = str(project)+'-'+str(data['run_accession'])+'-'+str(analysis_id))
def parse_samples(sample_sheet_path,connection,package,experiment_type):
    print ('Start Samples')
    sample_sheet_file = open(sample_sheet_path)
    column_names = sample_sheet_file.readline().strip('\n').split(',')
    to_add = []
    for line in sample_sheet_file:
        line = line.strip()
        project = line.split(',')[1]
        sample_name = line.split(',')[2]
        input_file_1_path = line.split(',')[3]
        input_file_2_path = line.split(',')[4]
        samples_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id)
        data = {'id':samples_id,'project':project,'sample_name':sample_name,'input_file_1_path':input_file_1_path,'input_file_2_path':input_file_2_path,
                'analysis_id':analysis_id,'experiment_type':experiment_type}
        if len(input_file_2_path)==0:
            del(data['input_file_2_path'])
            data['sequence_type'] = 'single'
        else:
            data['sequence_type'] = 'paired'
        to_add.append(data)
    #try:
    added_ids = add_multiple_rows(entity=package+'Samples',data=to_add,connection=connection)
    #except requests.exceptions.HTTPError as e:
    #    try:
    #        if 'Duplicate value' in e.response.json()['errors'][0]['message']:
    #            print(e.response.json()['errors'][0]['message'])
    #    except ValueError:
    #        raise
       
def parse_rnaseq_tools(sh_file_path,connection,package):
    '''filepath to .sh file used to run tool'''
    def time_from_log(logfile_text):
        '''.out file text with date stamps, returns runtime in seconds'''
        start_time_str = re.search('## \w+ (\w+ \d+ \d+:\d+:\d+ CEST \d{4}) [##|Start]',logfile_text).group(1)
        try:
            done_time_str = re.search('## \w+ (\w+ \d+ \d+:\d+:\d+ CEST \d{4}) ## Done',logfile_text).group(1)
        except:
            done_time_str = start_time_str
        start_time = datetime.datetime.strptime(start_time_str,'%b %d %H:%M:%S CEST %Y')
        done_time = datetime.datetime.strptime(done_time_str, '%b %d %H:%M:%S CEST %Y')
        delta_time = (done_time - start_time).total_seconds()
        return delta_time
    sh_id = None
    err_id = None
    out_id = None
    if not os.path.exists(os.path.dirname(sh_file_path)):
        raise OSError('Folder %s does not exist' % os.path.dirname(sh_file_path))
    sh_files = glob.glob(os.path.join(sh_file_path))
    if len(sh_files) == 0:
        raise ValueError ('For %s no files are found' % (sh_file_path))
    for sh_file in sh_files:
        sh_text = open(sh_file,'rb').read().decode("utf-8")
        split_sh = re.split('(## \S+ \S+ \d+:\d+:\d+ CEST \d+ ## \S+/slurm_script Started)',sh_text)
        if len(split_sh) > 2:
            sh_text = split_sh[-1] 

        tools = re.findall('module load (\S+)/(\S+)',sh_text)
        to_add = []
        for tool in tools:
            name = tool[0]
            version = tool[1]
            to_add.append({'id':name+'-'+version,'tool_name':name,'version':version})
        tool_ids = ','.join(add_multiple_rows(entity=package+'Tools',data=to_add, 
                                connection=connection))          
        basefile = os.path.splitext(sh_file)[0]
        err_text = open(basefile+'.err','rb').read().decode("utf-8") 
        out_text = open(basefile+'.out','rb').read().decode("utf-8") 
        try:
            runtime = str(int(time_from_log(out_text)))
        except AttributeError:
            warnings.warn(sh_file+' did not finish running (no done time found)')
            continue
        try:
            sample_name = re.search('sampleName="(.*?)"',sh_text).group(1)
        except AttributeError:
            sample_name = None
        try:
            internalId = re.search('internalId="(.*?)"',sh_text).group(1)
        except AttributeError:
            # easier to return None than only return internalId when it is found, None won't get added to the database either way
            internalId = None
        project = re.search('project="(.*?)"',sh_text).group(1)
        if not sh_id:
            sh_id = connection.add_file(file_path=basefile+'.sh', description='.sh script that was used to get the data', entity=package+'File',
                            data={'sample_file_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)+sh_file.replace(' ','_')})

        if not err_id:
            err_id = connection.add_file(file_path=basefile+'.err', description='file with messages printed to stderr by program', entity=package+'File',
                            data={'sample_file_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)+(basefile+'.err').replace(' ','_')})
        if not out_id:
            out_id = connection.add_file(file_path=basefile+'.out', description='file with messages printed to stdout by program',entity=package+'File',
                            data={'sample_file_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)+(basefile+'.out').replace(' ','_')})

        yield sh_text, err_text, out_text, runtime, sample_name, internalId, project, sh_id, err_id, out_id, tool_ids
def parse_depth_or_self(entity,depth_or_self_file,file_type,connection,package):
    '''parse depthRG or depthSM file'''
    depthIDs = []
    with open(depth_or_self_file) as verifyBamID_depth_or_self:
        verifyBamID_depth_or_self.readline()
        to_add = []
        count = 0
        for line in verifyBamID_depth_or_self:
            count += 1
            s_l = line.strip().split()
            if file_type == 'depth':
                data = {'rg':s_l[0],'depth':s_l[1],'number_of_SNPs':s_l[2],'perc_of_SNPs':s_l[3],'perc_cummulative':s_l[4]}
            elif file_type == 'self':
                data = {'seq_ID':s_l[0],'rg':s_l[1],'chip_id':s_l[2],'n_snps':s_l[3],'n_reads':s_l[4],'avg_dp':s_l[5],
                    'freemix':s_l[6],'freelk1':s_l[7],'freelk0':s_l[8],'free_ra':s_l[9],'chipmix':s_l[10],'chiplk1':s_l[11],'chiplk0':s_l[12],
                    'chip_rh':s_l[13],'chip_ra':s_l[14],'dpref':s_l[15],'rdphet':s_l[16],'rdpalt':s_l[17]}
            data_sanitized = {}
            for key, value in data.items():
                if value != 'NA' and value != 'NaN':
                    data_sanitized[key] = value
            to_add.append(data_sanitized)
        depthIDs = add_multiple_rows(package+entity, to_add, connection, False)
    return depthIDs
def parse_verifyBamID(runinfo_folder_QC,connection,package):
    '''finished'''
    print ('Start verifyBAMID')
    for sh_text, err_text, out_text, runtime,sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_QC,'VerifyBamID*.sh'),connection, package):
        verifyBamID_output_path = re.search('--out (\S+)',sh_text).group(1)
        with open(verifyBamID_output_path+'.log') as verifyBamID_log:
            verifyBamID_log_text = verifyBamID_log.read()
            groups = re.search('finding sample ID (\S+) from VCF file.*?'\
                              +'Finished reading (\d+) markers from VCF file.*?'\
                              +'Total of (\d+) informative markers passed.*?'\
                              +'Finished extracting (\d+) bases.*?'\
                              +'Avg Depth = (\d+\.\d+).*?',verifyBamID_log_text,re.DOTALL)
            number_of_markers = groups.group(2)
            number_of_informative_markers = groups.group(3)
            extracted_bases = groups.group(4) 
            avg_depth = groups.group(5)
            self_only_sample_ID = groups.group(1)
            self_only = bool(re.search('selfOnly option applied',out_text))
            to_add = []
            for skipped_marker_line in re.findall('Skipping.*?\n',verifyBamID_log_text):
                if len(skipped_marker_line.strip()) == 0:
                    continue
                groups = re.search('marker (\S+?):(\d+)',skipped_marker_line,re.DOTALL)
                chromosome = groups.group(1).replace('X','23')
                chromosome = chromosome.replace('Y','24')
                marker_type = 'no_autosomal' if 'no-autosomal' in skipped_marker_line else 'multiple_allele'
                position = groups.group(2)
                to_add.append({'chromosome':chromosome,'type':marker_type,'position':position})
            skipped_marker = ','.join(add_multiple_rows(entity=package+'Skipped_marker',data=to_add,
                                                        connection=connection))
            pattern = '(Comparing with individual \S+.. Optimal fIBD = \d+\.\d+, LLK0 = \d+\.\d+, LLK1 = \d+\.\d+ for readgroup \d+\n'\
                     +'Best Matching Individual is NA with IBD = \d+.\d+\n'\
                     +'Self Individual is NA with IBD = \d+.\d+)'
            verifyBamID_individual = ''
            to_add = []
            for individual in re.findall(pattern, verifyBamID_log_text):
                groups = re.search('Optimal fIBD = (\d+\.\d+).*?'\
                                  +'LLK0 = (\d+\.\d+).*?'\
                                  +'LLK1 = (\d+\.\d+).*?'\
                                  +'for readgroup (\S+)\n.*?'\
                                  +'Best Matching Individual is (\S+) with IBD.*?',
                                  +'Best Matching Individual is \S+ with IBD = (\d+.\d+)\n.*?'\
                                  +'Self Individual is (\S+) with IBD.*?'\
                                  +'Self Individual is \S+ with IBD = (\d+.\d+)',individual,re.DOTALL)
                data = {'optimal_fIBD':groups.group(1),'llk0':groups.group(2),'llk1':groups.group(3),'readgroup':groups.group(4),'best_matching_individual':groups.group(5),
                        'best_match_individual_IBD':groups.group(6),'self_individual':groups.group(7),'self_individual_IBD':groups.group(8)}
                to_add.append(data)
            verifyBamID_individual = ','.join(add_multiple_rows(package+'VerifyBamID_individual', to_add,connection,False))
        depthRG = ','.join(parse_depth_or_self('DepthRG',verifyBamID_output_path+'.depthRG','depth', connection,package))
        depthSM = ','.join(parse_depth_or_self('DepthSM',verifyBamID_output_path+'.depthSM','depth', connection,package))
        selfRG = ','.join(parse_depth_or_self('SelfRG',verifyBamID_output_path+'.selfRG','self', connection,package))
        selfSM = ','.join(parse_depth_or_self('SelfSM',verifyBamID_output_path+'.selfSM','self', connection,package))
        multiple_allele_skip = len(re.findall('with multiple alternative alleles', err_text))
        no_autosomal_skip = len(re.findall('Skipping no-autosomal marker',err_text))
        data = {'multiple_allele_skip':multiple_allele_skip,'no_autosomal_skip':no_autosomal_skip,'sh_script':sh_id, 'depthRG':depthRG,
                'out_file':out_id,'err_file':err_id,'runtime':runtime, 'depthSM':depthSM,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'number_of_markers':number_of_markers,'number_of_informative_markers':number_of_informative_markers,'verifyBamID_individual':verifyBamID_individual,
                'selfRG':selfRG,'selfSM':selfSM,'extracted_bases':extracted_bases,'avg_depth':avg_depth,'skipped_marker':skipped_marker,
                'self_only':self_only,'self_only_sample_ID':self_only_sample_ID, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}        
        added_id = connection.add(package+'VerifyBamID', data)[0]
        # if this sample already had a verifyBamId, find it, append the new one, and update sample row
        verifybamid_data = connection.get(package+'VerifyBamID', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(verifybamid_data) >0 and len(verifybamid_data[0]['id']) > 0:
            added_id = verifybamid_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'verifyBamID':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_hisat(runinfo_folder_QC,connection,package):
    '''finished'''
    print ('Start Hisat')
    for hisat_sh_text, hisat_err_text, hisat_out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_QC,'HisatAlignment*.sh'),connection, package):
        groups = re.search('(\d+) reads; of these:.*?'\
                          +'(\d+) \((\d+\.\d+)%\) were paired; of these:.*?'\
                          +'(\d+) \((\d+\.\d+)%\) aligned concordantly 0 times.*?'\
                          +'(\d+) \((\d+\.\d+)%\) aligned concordantly exactly 1 time.*?'\
                          +'(\d+) \((\d+\.\d+)%\) aligned concordantly \>1 times.*?'\
                          +'(\d+) \((\d+\.\d+)%\) aligned discordantly 1 time.*?'\
                          +'(\d+) pairs aligned 0 times concordantly or discordantly; of these:.*?'
                          +'(\d+) mates make up the pairs; of these:.*?'\
                          +'(\d+) \((\d+.\d+)%\) aligned 0 times.*?'\
                          +'(\d+) \((\d+\.\d+)%\) aligned exactly 1 time.*?'\
                          +'(\d+) \((\d+\.\d+)%\) aligned \>1 times.*?'\
                          +'(\d+\.\d+)% overall alignment rate', hisat_err_text,re.DOTALL)
        data = {'number_of_reads':groups.group(1),'number_of_paired_reads':groups.group(2),'number_of_paired_reads_p':groups.group(3),'align_pairs_conc_0_times':groups.group(4),'align_pairs_conc_0_times_p':groups.group(5),'align_pairs_conc_1_time':groups.group(6),'align_pairs_conc_1_time_p':groups.group(7),
                'align_pairs_conc_2plus_time':groups.group(8),'align_pairs_conc_2plus_time_p':groups.group(9),'aligned_pairs_disconc':groups.group(10),'aligned_pairs_disconc_p':groups.group(11),'pairs_not_aligned':groups.group(12),
                'number_mates':groups.group(13),'mates_aligned_0_times':groups.group(14),'mates_aligned_0_times_p':groups.group(15),'mates_aligned_1_time':groups.group(16),'mates_aligned_1_time_p':groups.group(17),
                'mates_aligned_multiple':groups.group(18),'mates_aligned_multiple_p':groups.group(19),'overall_alignment_rate':groups.group(20),'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'sh_script':sh_id,'out_file':out_id,'err_file':err_id,'runtime':runtime, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'Hisat', data)[0]
        hisat_data = connection.get(package+'Hisat', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(hisat_data) >0 and len(hisat_data[0]['id']) > 0:
            added_id = hisat_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'hisat':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_variantEval(runinfo_folder_QC,connection,package):
    '''finished'''
    print ('Start variantEval')
    for sh_text, err_text, text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_QC,'VariantEval_*.sh'), connection,package):
        table_name = None
        table_name_previous = None
        variant_eval_ids = {'indel_summary':'','comp_overlap':'','variant_class_count':'',
                            'indel_length_histogram':'','multiallelic_summary':'','ti_tv_variant_evaluator':'',
                            'validation_report':'','variant_summary':''}
        variant_eval_file_path = re.search('-o (\S+.grp)',sh_text).group(1)
        # remove below line when on cluster
        #variant_eval_file_path = '/Users/Niek/UMCG/test/data/ATACseq/project/variantEval/'+variant_eval_file_path.split('/')[-1]
        with open(variant_eval_file_path) as variant_eval_file:
            to_add = []
            for line in variant_eval_file:
                line = line.strip()
                if len(line) == 0:
                    continue
                s_l = line.split()
                table_name = s_l[0]
                if table_name == table_name_previous:
                    try:
                        if table_name == 'IndelSummary':
                            data = {'comp_rod':s_l[1],'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'n_SNPs':s_l[5],'n_singletons':s_l[6],'n_indels':s_l[7],'n_singleton_indels':s_l[8],
                                    'n_indels_match_gold_standard':s_l[9],'gold_standard_match_rate':s_l[10],'n_multiallelic_indel_sites':s_l[11],'perc_3_or_more_alleles':s_l[12],
                                    'snp_to_indel_ratio':s_l[13],'snp_indel_ratio_singletons':s_l[14], 'n_novel_indels':s_l[15], 'indel_novelty_rate':s_l[16], 'n_insertions':s_l[17], 'n_dels':s_l[18],
                                    'insertion_to_del_ratio':s_l[19],'n_large_dels':s_l[20],'n_large_insertions':s_l[21],'insert_del_ratio_large_indels':s_l[22], 'n_coding_indels_frameshifting':s_l[23],
                                    'n_coding_indels_in_frame':s_l[24],'frameshift_rate_coding_indels':s_l[25],'snp_het_to_hom_ratio':s_l[26], 'indel_het_to_hom_ratio':s_l[27], 'ratio_1_2_to_3_bp_insertions':s_l[28],
                                    'ratio_1_2_to_3_bp_dels':s_l[29]}
                        elif table_name == 'CompOverlap':
                            data = {'comp_rod':s_l[1], 'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'n_eval_variants':s_l[5],
                                    'novel_sites':s_l[6],'n_variants_at_comp':s_l[7],'comp_rate':s_l[8],'n_concant':s_l[9],'concant_rate':s_l[10]}
                        elif table_name == 'CountVariants':
                            data = {'comp_rod':s_l[1],'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'n_processed_loci':s_l[5],'n_called_loci':s_l[6],'n_ref_loci':s_l[7],'n_variant_loci':s_l[8],
                                    'variant_rate':s_l[9],'variant_rate_per_bp':s_l[10],'n_SNPs':s_l[11],'n_MNPs':s_l[12],'n_insertions':s_l[13],'n_dels':s_l[14],'n_complex':s_l[15],'n_symbolic':s_l[16],
                                    'n_mixed':s_l[17],'n_no_calls':s_l[18],'n_hets':s_l[19],'n_hom_ref':s_l[20],'n_hom_var':s_l[21],'n_singletons':s_l[22],'n_hom_derived':s_l[23],'heterozygosity':s_l[24],
                                    'heterozygosity_per_bp':s_l[25],'het_hom_ratio':s_l[26],'indel_rate':s_l[27],'indel_rate_per_bp':s_l[28],'insertion_del_ratio':s_l[29]}
                            table_name = 'Variant_class_count'
                        elif table_name == 'IndelLengthHistogram':
                            data = {'comp_rod':s_l[1],'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'length':s_l[5],'freq':s_l[6]}
                        elif table_name == 'MultiallelicSummary':
                            data = {'comp_rod':s_l[1],'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'n_processed_loci':s_l[5],'n_SNPs':s_l[6],'n_multi_SNPs':s_l[7],'processed_multi_snp_ratio':s_l[8],
                                    'variant_multi_snp_ratio':s_l[9],'n_indels':s_l[10],'n_multi_indels':s_l[11],'processed_multi_indel_ratio':s_l[12],'variant_multi_indel_ratio':s_l[13],'n_ti':s_l[14],'n_tv':s_l[15],
                                    'ti_tv_ratio':s_l[16],'known_SNPs_partial':s_l[17],'known_SNPs_complete':s_l[18],'snp_novelty_rate':s_l[19]}
                        elif table_name == 'TiTvVariantEvaluator':
                            data = {'comp_rod':s_l[1],'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'n_ti':s_l[5],'n_tv':s_l[6],'ti_tv_ratio':s_l[7],'n_ti_in_comp':s_l[8],'n_tv_in_comp':s_l[9],
                                    'ti_tv_ratio_standard':s_l[10],'n_ti_derived':s_l[11],'n_tv_derived':s_l[12],'ti_tv_derived_ratio':s_l[13]}
                        elif table_name == 'ValidationReport':
                            data = {'comp_rod':s_l[1],'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'n_comp':s_l[5],'tp':s_l[6],'fp':s_l[7],'fn':s_l[8],'tn':s_l[9],'sensitivity':s_l[10],
                                    'specificity':s_l[11],'ppv':s_l[12],'fdr':s_l[13],'comp_mono_eval_no_call':s_l[14],'comp_mono_eval_filtered':s_l[15],'comp_mono_eval_mono':s_l[16],'comp_mono_eval_poly':s_l[17],
                                    'comp_poly_eval_no_call':s_l[18],'comp_poly_eval_filtered':s_l[19],'comp_poly_eval_mono':s_l[20],'comp_poly_eval_poly':s_l[21],'comp_filtered':s_l[22],'n_different_allele_sites':s_l[23]}
                        elif table_name == 'VariantSummary':
                            data = {'comp_rod':s_l[1],'eval_rod':s_l[2],'jexl_expression':s_l[3],'novelty':s_l[4],'n_samples':s_l[5],'n_processed_loci':s_l[6],'n_SNPs':s_l[7],'ti_tv_ratio':s_l[8],'snp_novelty_rate':s_l[9],
                                    'n_SNPs_per_sample':s_l[10],'ti_tv_ratio_per_sample':s_l[11],'snp_DP_per_sample':int(float(s_l[12])),'n_indels':s_l[13],'indel_novelty_rate':s_l[14],'n_indels_per_sample':s_l[15],'indel_DP_per_sample':int(float(s_l[16])),
                                    'n_SVs':s_l[17],'sv_novelty_rate':s_l[18],'n_SVs_per_sample':s_l[19]}
                    except:
                        print (line)
                        print (s_l)
                        print((len(s_l)))
                        raise
                    data_sanitized = {}
                    for key, value in data.items():
                        if value != 'NA' and value != 'NaN':
                            data_sanitized[key] = value
                    to_add.append(data_sanitized)
                    # table name to entity name example: VariantSummary to Variant_summary
                    entity = re.sub( '(?<!^)(?=[A-Z])', '_', table_name ).lower() 
                table_name_previous = table_name
            variant_eval_ids[entity.lower()] += add_multiple_rows(entity=package+entity,data=to_add,
                                                                                connection=connection)[0]

        data = {'indel_summary':variant_eval_ids['indel_summary'], 'variant_comp_overlap':variant_eval_ids['comp_overlap'],'indel_length_histogram':variant_eval_ids['indel_length_histogram'],
                'multiallelic_summary':variant_eval_ids['multiallelic_summary'],'ti_tv_variant_evaluator':variant_eval_ids['ti_tv_variant_evaluator'],'validation_report':variant_eval_ids['validation_report'],
                'variant_summary':variant_eval_ids['variant_summary'], 'out_file':out_id,'err_file':err_id, 'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'runtime':runtime,'variant_class_count':variant_eval_ids['variant_class_count'],'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        for key, value in data.items():
            data[key] = value.rstrip(',')
        added_id = connection.add(package+'VariantEval',data)[0]
        variantEval_data = connection.get(package+'hisat', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(variantEval_data) >0 and len(variantEval_data[0]['id']) > 0:
            added_id = variantEval_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'variantEval':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_md5sums(connection, package):
    '''TODO: LINK TO SAMPLES!!!!'''
    print ('start md5sums')
    project_folder = configSectionMap("paths")['project_folder']
    x = 0
    for root, dirs, files in os.walk(project_folder):
        for name in files:
            if name.endswith('.md5'):
                x += 1
                for line in open(root+'/'+name):
                    data = {'md5sum':line.split('  ')[0], 'file_name':line.split('  ')[1]}
                    connection.add(package+'Md5sums',data)
    if x == 0:
        print('No .md5 files found in %s or its subdirectories' % project_folder)
def parse_variantCaller(variant_caller, runinfo_folder_QC,connection,package):
    print('Start variantCaller')
    pretext = ''
    if variant_caller != 'UnifiedGenotyper' and variant_caller != 'HaplotypeCaller' and variant_caller != 'GenotypeGvcf':
        raise AttributeError('variant_caller not UnifiedGenotyper, HaplotypeCaller or GenotypeGvcf')
    if variant_caller == 'UnifiedGenotyper' or variant_caller == 'HaplotypeCaller':
        pretext = 'Gatk'
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_QC,pretext+variant_caller+'*.sh'),connection,package):
        if variant_caller != 'GenotypeGvcf':
            total_reads = re.search('out of approximately (\d+) total reads',err_text).group(1)
            reads_filtered_out = re.search('(\d+) reads were filtered out',err_text).group(1)
            reads_filtered_out_perc = re.search('out of approximately \d+ total reads \((\d+.\d+)%\)',err_text).group(1)
            filtered = {}
            filters = ['BadCigarFilter','BadMateFilter','DuplicateReadFilter','FailsVendorQualityCheckFilter','MalformedReadFilter',
                       'MappingQualityUnavailableFilter','NotPrimaryAlignmentFilter','ReassignMappingQualityFilter','UnmappedReadFilter']
            if variant_caller == 'HaplotypeCaller':
                filters.append('HCMappingQualityFilter')
                
            for reads_filter in filters:
                search_result = re.search('(\d+) reads \((\d+.\d+)% of total\) failing '+reads_filter,err_text)
                filtered[reads_filter] = {'p':search_result.group(2),'r':search_result.group(1)}
        vcf_file = re.search('-o (\S+.vcf)',sh_text).group(1)
        sample_names = None
        with open(vcf_file) as vcf:
            vcf_meta_ids = {'filters':[],'info':[],'formats':[],'contigs':[],'snps':[],'alts':[],'gvcf_block':[]}
            reference = None
            for line in vcf:
                line = line.strip()
                if line.startswith('#'):
                    if '##contig' in line:
                        search_filter = re.search('ID=(\S+),length=(\d+),assembly=(\S+)>',line)
                        vcf_meta_ids['contigs'] += [{'meta_id':search_filter.group(1),'length':search_filter.group(2),'assembly':search_filter.group(3)}]
                    elif '##INFO' in line:
                        search_filter = re.search('ID=(\S+),Number=(\S),Type=(\w+),Description="(.*?)">',line)
                        vcf_meta_ids['info'] += [{'meta_id':search_filter.group(1),'number':search_filter.group(2),'type':search_filter.group(3).lower(),'description':search_filter.group(4)}]
                    elif '##FORMAT' in line:
                        search_filter = re.search('ID=(\S+),Number=(\S),Type=(\w+),Description="(.*?)">',line)
                        vcf_meta_ids['formats'] += [{'meta_id':search_filter.group(1),'number':search_filter.group(2),'type':search_filter.group(3).lower(),'description':search_filter.group(4)}]
                    elif '##FILTER' in line:
                        search_filter = re.search('ID=(\S+),Description="(.*?)">',line)
                        vcf_meta_ids['filters'] += [{'meta_id':search_filter.group(1),'description':search_filter.group(2)}]
                    elif '##ALT' in line:
                        search_filter = re.search('ID=(\S+),Description="(.*?)">',line)
                        vcf_meta_ids['alts'] += [{'meta_id':search_filter.group(1),'description':search_filter.group(2)}]
                    elif '##GVCFBlock' in line:
                        search_filter = re.search('GVCFBlock(\d+-\d+)=(\d+)\S+maxGQ=(\d+)',line)
                        vcf_meta_ids['gvcf_block'] += [{'gvcf_block_min':search_filter.group(1),'gvcf_block_max':search_filter.group(2), 'min_gq':search_filter.group(3),'max_gq':search_filter.group(4)}]
                    elif 'fileformat=' in line:
                        fileformat = line.split('fileformat=')[1]
                    elif 'reference=' in line:
                        reference = line.split('reference=')[1]
                    elif '#CHROM' in line:
                        sample_names = line.split('\t')[9:]
                    else:
                        break
            #### Too much data in vcf file, cant use commented data
            #        s_l = line.split('\t')
            #        snp_per_sample = s_l[9:]
            #        snp_per_sample_id = ''
            #        for sample_name, snp_info in zip(sample_names, snp_per_sample):
            #            snp_per_sample_id += connection.add(package+'public_rnaseq_Snp_per_sample',data = {'sample_name':sample_name,'snp_info':snp_info})+','
            #        snp_per_sample_id = snp_per_sample_id.rstrip(',')
            #        vcf_meta_ids['snps'] += connection.add(package+'Snp_info',{'chrom':s_l[0].replace('X','23').replace('Y','24'),'pos':s_l[1],'id':s_l[2],'ref':s_l[3],'alt':s_l[4],'qual':s_l[5],'filter':s_l[6],'info':s_l[7],'format':s_l[8],'snp_per_sample':snp_per_sample_id})+','
            if vcf_meta_ids['contigs']:
                contig_ids = ','.join(add_multiple_rows(entity=package+'Contigs',data=vcf_meta_ids['contigs'], connection=connection))
            else:
                contig_ids = None
            if vcf_meta_ids['info']:
                info_id = add_multiple_rows(entity=package+'Info',data=vcf_meta_ids['info'], connection=connection)
            else:
                info_id = None
            if vcf_meta_ids['formats']:
                format_ids = ','.join(add_multiple_rows(entity=package+'Formats',data=vcf_meta_ids['formats'], connection=connection))
            else:
                format_ids = None
            if vcf_meta_ids['filters']:
                filter_ids = ','.join(add_multiple_rows(entity=package+'Filters',data=vcf_meta_ids['filters'], connection=connection))
            else:
                filter_ids = None
            if vcf_meta_ids['alts']:
                alt_ids = ','.join(add_multiple_rows(entity=package+'Alts',data=vcf_meta_ids['alts'], connection=connection))
            else:
                alt_ids = None
            if vcf_meta_ids['gvcf_block']:
                gvcf_blck_ids = ','.join(add_multiple_rows(entity=package+'Gvcf_block',data=vcf_meta_ids['gvcf_block'], connection=connection))
            else:
                gvcf_blck_ids = None
            data = {'fileformat':fileformat,'contigs':contig_ids,'filters':filter_ids,'formats':format_ids,'reference':reference,
                    'alts':alt_ids,'gvcf_block':gvcf_blck_ids,'info':info_id}
            added_id = connection.add(package+'Vcf',data)
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id),'type':variant_caller} # 'vcf':added_id,
        if variant_caller != 'GenotypeGvcf':
            extra_data = {'unmapped_read_filter_perc':filtered['UnmappedReadFilter']['p'],'unmapped_read_filter':filtered['UnmappedReadFilter']['r'],'total_reads':total_reads,'reads_filtered_out_perc':reads_filtered_out_perc,'reads_filtered_out':reads_filtered_out,
                          'not_primary_alignment_filter':filtered['NotPrimaryAlignmentFilter']['r'],'not_primary_align_filter_perc':filtered['NotPrimaryAlignmentFilter']['p'],'map_qual_unavailable_filter':filtered['MappingQualityUnavailableFilter']['r'],'map_qual_unavail_filter_perc':filtered['MappingQualityUnavailableFilter']['p'],'malformed_read_filter_perc':filtered['MalformedReadFilter']['p'],
                          'malformed_read_filter':filtered['MalformedReadFilter']['r'],'fails_vendor_qual_filter_perc':filtered['FailsVendorQualityCheckFilter']['p'],'fails_vendor_qual_filter':filtered['FailsVendorQualityCheckFilter']['r'],'duplicate_read_filter_perc':filtered['DuplicateReadFilter']['p'],'duplicate_read_filter':filtered['DuplicateReadFilter']['r'],
                          'bad_mate_filter_perc':filtered['BadMateFilter']['p'],'bad_mate_filter':filtered['BadMateFilter']['r'],'bad_cigar_filtered_perc':filtered['BadCigarFilter']['p'],'bad_cigar_filtered':filtered['BadCigarFilter']['r']}
            data.update(extra_data)
        if variant_caller == 'HaplotypeCaller':
            data.update({'hc_mapping_quality_filter':filtered['HCMappingQualityFilter']['r'],'hc_mapping_quality_perc':filtered['HCMappingQualityFilter']['p']})
        if variant_caller == 'UnifiedGenotyper':
            data.update({'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId})
        added_id = connection.add(package+'VariantCaller', data)[0]
        if variant_caller == 'UnifiedGenotyper':
            unifiedGenotyper_data = connection.get(package+'VariantCaller', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
            if len(unifiedGenotyper_data) >0 and len(unifiedGenotyper_data[0]['id']) > 0:
                added_id = unifiedGenotyper_data[0]['id']+','+added_id
        if sample_names:
            for sample_name in sample_names:
                connection.update_entity_rows(package+'Samples', data={'variantCaller':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
        else:
            connection.update_entity_rows(package+'Samples', data={'variantCaller':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_fastqc(runinfo_folder_QC,connection,package):
    '''todo: readlength can be 50 instead of 30-50, .rstrip to ','.join()'''
    print('Start fastqc')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_QC,'FastQC*.sh'),connection, package):
        fastqc_outdir = re.search('--outdir (\S+)',sh_text).group(1)
        fastqc_output_base = re.search('--noextract (\S+)',sh_text).group(1).rstrip('.fastq.gz').split(os.sep)[-1]+'_fastqc'
        with open(fastqc_outdir+os.sep+fastqc_output_base+'.html') as fastqc_html_file:
            fastqc_html = fastqc_html_file.read()
            images = re.findall('img src="data:image/png;base64,(\S+)" alt="(\S+)"/>\<a href="\S+"\>(.*?)\<',fastqc_html)
            image_data = {}
            for image in images:
                image_data[image[2]] = [image[0],image[1].strip('[').strip(']').lower()]
        with zipfile.ZipFile(fastqc_outdir+os.sep+fastqc_output_base+'.zip') as archive:
            graphs = {}
            names = []
            for name in archive.namelist():
                if 'Images' in name and '.png' in name:
                    description = os.path.basename(os.path.splitext(name)[0])
                    img_data = archive.open(name)
                    names.append(description)
                    graphs[description] = connection.add(name, description,package+'File', io_stream = img_data,
                                                              extra_data={'sample_file_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)+name.replace(' ','_')})
                if 'fastqc_data.txt' in name:
                    fastqc_data = archive.read(name)
        fastqc_groups = re.search('Filename\s+(\S+)\n'\
                                 +'File type\s+(.+)\n'\
                                 +'Encoding\s+(.+)\n'\
                                 +'Total Sequences\s+(\S+)\n'\
                                 +'Sequences flagged as poor quality\s+(\d+)\n'\
                                 +'Sequence length\s+(\d+-\d+)\n'\
                                 +'%GC\s+(\d+)\n.*?'\
                                 +'Per base sequence quality.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Per tile sequence quality.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Per sequence quality scores.*?\n(.*)\n\>\>END_MODULE.*?'\
                                 +'Per base sequence content.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Per sequence GC content.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Per base N content.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Sequence Length Distribution.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Sequence Duplication Levels.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Overrepresented sequences.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Adapter Content.*?\n(.*?)\n\>\>END_MODULE.*?'\
                                 +'Kmer Content.*?\n(.*?)\n\>\>END_MODULE', fastqc_data.decode('utf-8'),re.DOTALL)
        
        fastqc_ids = {'pbsq':'','ptsq':'','psqs':'','pbsc':'','psgc':'','pbnc':'','sld':'','sdl':'','os':'','ac':'','kc':''}
        to_add = []
        for line in fastqc_groups.group(8).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'base_letter':s_l[0], 'mean':s_l[1], 'median':s_l[2], 'lower_quartile':s_l[3], 'upper_quartile':s_l[4],'percentile_10th':s_l[5], 'percentile_90th':s_l[6]})
        fastqc_ids['pbsq'] = ','.join(add_multiple_rows(package+'Per_base_sequence_quality',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(9).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'tile':s_l[0], 'base_letter':s_l[1], 'mean':s_l[2]})
        fastqc_ids['ptsq'] = ','.join(add_multiple_rows(package+'Per_tile_sequence_quality',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(10).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'quality':s_l[0], 'count':s_l[1]})
        fastqc_ids['psqs'] = ','.join(add_multiple_rows(package+'Per_sequence_quality_scores',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(11).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'base_letter':s_l[0], 'a':s_l[1], 'g':s_l[2], 'c':s_l[3], 't':s_l[4]})
        fastqc_ids['pbsc'] = ','.join(add_multiple_rows(package+'Per_base_sequence_content',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(12).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'gc_content':s_l[0], 'count':s_l[1]})
        fastqc_ids['psgc'] = ','.join(add_multiple_rows(package+'Per_sequence_GC_content',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(13).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'base_letter':s_l[0], 'n_count':s_l[1]})
        fastqc_ids['pbnc'] = ','.join(add_multiple_rows(package+'Per_base_N_content',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(14).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'length':s_l[0], 'count':s_l[1]})
        fastqc_ids['sld'] = ','.join(add_multiple_rows(package+'Sequence_length_distribution',to_add,connection,True))
        total_deduplicated_perc = fastqc_groups.group(15).split('\n')[0].split('\t')[1]            
        to_add = []
        for line in fastqc_groups.group(15).split('\n')[2:]:
            s_l = line.split('\t')
            to_add.append({'duplication_level':s_l[0],'perc_of_deduplicated':s_l[1],'perc_of_total':s_l[2]})
        fastqc_ids['sdl'] = ','.join(add_multiple_rows(package+'Sequence_duplication_levels',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(16).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'sequence':s_l[0],'count':s_l[1],'percentage':s_l[2],'possible_source':s_l[3]})
        fastqc_ids['os'] = ','.join(add_multiple_rows(package+'Overrepresented_sequences',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(17).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'position':s_l[0], 'illumina_universal_adapter':s_l[1], 'illumina_small_rna_adapter':s_l[2], 'nextera_transposase_seq':s_l[3], 'solid_small_rna_adapter':s_l[4]})
        fastqc_ids['ac'] = ','.join(add_multiple_rows(package+'Adapter_content',to_add,connection,True))
        to_add = []
        for line in fastqc_groups.group(18).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'seq':s_l[0], 'count':s_l[1], 'p_value':s_l[2], 'obs_exp_max':s_l[3], 'max_obs_exp_position':s_l[4]})
        fastqc_ids['kc'] = ','.join(add_multiple_rows(package+'Kmer_content',to_add,connection,True))
        if '-' in fastqc_groups.group(6):
            seq_length_min = fastqc_groups.group(6).split('-')[0]
        else:
            seq_length_min = fastqc_groups.group(6)
        if '-' in fastqc_groups.group(6):
            seq_length_max = fastqc_groups.group(6).split('-')[1]
        else:
            seq_length_max = fastqc_groups.group(6)
        data = {'adapter_content_graph':graphs['adapter_content'],'kmer_content_graph':graphs['kmer_profiles'],'per_base_n_content_graph':graphs['per_base_n_content'],
                'per_base_seq_content_graph':graphs['per_base_sequence_content'],'per_base_seq_qual_graph':graphs['per_base_quality'],'per_seq_GC_content_graph':graphs['per_sequence_gc_content'],
                'per_seq_qual_scores_graph':graphs['per_sequence_quality'],'per_tile_seq_qual_graph':graphs['per_tile_quality'],'seq_duplication_levels_graph':graphs['duplication_levels'],
                'seq_length_distribution_graph':graphs['sequence_length_distribution'],
                'total_deduplicated_perc':total_deduplicated_perc,'adapter_content':fastqc_ids['ac'],'adapter_content_check':image_data['Adapter Content'][1],'basic_statistics_check':image_data['Basic Statistics'][1],'gc_perc':fastqc_groups.group(7),
                'kmer_content':fastqc_ids['kc'],'kmer_content_check':image_data['Kmer Content'][1],'overrepresented_seqs':fastqc_ids['os'],
                'overrepresented_seqs_check':image_data['Overrepresented sequences'][1],'per_base_N_content':fastqc_ids['pbnc'],'per_base_n_content_check':image_data['Per base N content'][1],'per_base_seq_content':fastqc_ids['pbsc'],
                'per_base_seq_content_check':image_data['Per base sequence content'][1],'per_base_seq_qual':fastqc_ids['pbsq'],
                'per_base_seq_qual_check':image_data['Per base sequence quality'][1],'per_seq_GC_content':fastqc_ids['psgc'],'per_seq_gc_content_check':image_data['Per sequence GC content'][1],
                'per_seq_qual_scores':fastqc_ids['psqs'],'per_seq_qual_scores_check':image_data['Per sequence quality scores'][1],'per_tile_seq_qual':fastqc_ids['ptsq'],
                'per_tile_seq_qual_check':image_data['Per tile sequence quality'][1],'seq_duplication_levels':fastqc_ids['sdl'],
                'seq_duplication_levels_check':image_data['Sequence Duplication Levels'][1],'seq_length_min':seq_length_min,'seq_length_distribution':fastqc_ids['sld'],
                'seq_length_distribution_check':image_data['Sequence Length Distribution'][1],'seqs_flagged_as_poor':fastqc_groups.group(5),'total_seqs':fastqc_groups.group(4),'seq_length_max':seq_length_max,'file_name':fastqc_groups.group(2),'file_type':fastqc_groups.group(2),'encoding':fastqc_groups.group(3),
                'err_file':err_id,'out_file':out_id,'runtime':runtime,'sh_script':sh_id,'internalId':internalId,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name)}
        added_id = connection.add(package+'FastQC', data)[0]
        fastqc_data = connection.get(package+'fastqc', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(fastqc_data) >0 and len(fastqc_data[0]['id']) > 0:
            added_id = fastqc_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'fastqc':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_markDuplicates(runinfo_folder_genotypeCalling,connection,package):
    print ('Start markDuplicates')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'MarkDuplicates*.sh'),connection, package):
        duplicates = re.search('Marking (\d+) records as duplicates',err_text).group(1)
        optical_duplicate_clusters = re.search('Found (\d+) optical duplicate clusters',err_text).group(1)
        mark_duplicates_log = re.search('OUTPUT=(.*?)\.bam',sh_text).group(1)+'.metrics.log'
        # remove below line when on cluster
        #mark_duplicates_log ='/Users/Niek/UMCG/test/data/ATACseq/project/markDuplicates/'+mark_duplicates_log.split('/')[-1]
        duplication_metrics = ''
        markduplicates_histogram = ''
        with open(mark_duplicates_log) as log_file:
            log = log_file.read()
            groups = re.search('METRICS CLASS(.*?)## HISTOGRAM(.*?)',log,re.DOTALL)
            to_add = []
            for line in groups.group(1).split('\n'):
                if not line.startswith('LIBRARY') and not line.startswith('##') and len(line.split('\t')[0].strip()) > 0:
                    s_l = line.split('\t')
                    to_add.append({'library':s_l[0],'unpaired_reads_examined':s_l[1],'reads_pairs_examined':s_l[2],'unmapped_reads':s_l[3],'unpaired_read_duplicates':s_l[4],
                                   'read_pair_duplicates':s_l[5],'read_pair_optical_duplicates':s_l[6],'percent_duplication':s_l[7],'estimated_library_size':s_l[8]})
            duplication_metrics = ','.join(add_multiple_rows(package+'Duplication_metrics', to_add,connection,True))
            to_add = []
            for line in groups.group(2).split('\n'):
                if not line.startswith('##') and not line.startswith('BIN') and len(line.strip()) > 0:
                    to_add.append({'bin':line.split('\t')[0],'value':line.split('\t')[1]})
            markduplicates_histogram = ','.join(add_multiple_rows(package+'MarkDuplicates_histogram',to_add,connection,True))

        data = {'duplicates':duplicates,'optical_duplicate_clusters':optical_duplicate_clusters,
                'duplication_metrics':duplication_metrics.rstrip(','), 'markduplicates_histogram':markduplicates_histogram.rstrip(','),
                'err_file':err_id,'out_file':out_id,'runtime':runtime,'sh_script':sh_id}
        added_id = connection.add(package+'MarkDuplicates', data)[0]
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'markDuplicates':added_id})
def parse_bqsr(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print ('start BQSR')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'BQSR*.sh'),connection, package):
        bqsr_output_folder = re.search('module list.*?mkdir -p (\S+)',sh_text,re.DOTALL).group(1)
        # remove below line when on cluster
        #bqsr_output_folder ='/Users/Niek/UMCG/test/data/ATACseq/project/baseQualityScoreRecalibration/'
        before_after_ids = {}
        for before_after in ['before','after']:
            with open(bqsr_output_folder+sample_name+'.'+before_after+'.grp') as grp:
                arguments_table = False
                quantized_table = False
                recal_table_0 = False
                recal_table_1 = False
                recal_table_2 = False
                arguments = ''
                quantized = ''
                recal_0 = ''
                recal_1 = ''
                recal_2 = ''
                skip_line = False
                to_add_dict = {'Arguments':[],'Quantized':[],'Recal_table_0':[],'Recal_table_1':[],'Recal_table_2':[]}
                for line in grp:
                    if skip_line:
                        skip_line = False
                        continue
                    line = line.strip()
                    if len(line) == 0:
                        arguments_table = False
                        quantized_table = False
                        recal_table_0 = False
                        recal_table_1 = False
                        recal_table_2 = False
                    s_l = line.split()
                    if arguments_table:
                        to_add_dict['Arguments'].append({'argument':s_l[0],'value':s_l[1]})
                    elif quantized_table:
                        to_add_dict['Quantized'].append({'quality_score':s_l[0],'count':s_l[1],'quantized_score':s_l[2]})
                    elif recal_table_0:
                        data = {'event_type':s_l[1],'empirical_quality':s_l[2],'estimated_q_reported':s_l[3],
                                'errors':s_l[4],'observations':s_l[5]}
                        to_add_dict['Recal_table_0'].append(data)
                    elif recal_table_1:
                        data = {'quality_score':s_l[1],'event_type':s_l[2],'empirical_quality':s_l[3],
                                'observations':s_l[4],'errors':s_l[5]}
                        to_add_dict['Recal_table_1'].append(data)
                    elif recal_table_2:
                        data = {'quality_score':s_l[1],'covariate_value':s_l[2],'covariate_name':s_l[3],
                                'event_type':s_l[4],'empirical_quality':s_l[5],'observations':s_l[6],'errors':s_l[7]}
                        to_add_dict['Recal_table_2'].append(data)
                        #recal_table_2 = False
                    elif 'GATKTable:Arguments' in line:
                        arguments_table = True
                        skip_line = True
                    elif 'GATKTable:Quantized' in line:
                        quantized_table = True
                        skip_line = True
                    elif 'GATKTable:RecalTable0' in line:
                        recal_table_0 = True
                        skip_line = True
                    elif 'GATKTable:RecalTable1' in line:
                        recal_table_1 = True
                        skip_line = True
                    elif 'GATKTable:RecalTable2' in line:
                        recal_table_2 = True
                        skip_line = True
            mref_ids = {}
            for key in to_add_dict:
                mref_ids[key] = ','.join(add_multiple_rows(package+key, to_add_dict[key],connection,True))
            data = {'argument':mref_ids['Arguments'],'quantized':mref_ids['Quantized'],'recal_table_1':mref_ids['Recal_table_1'],'recal_table_2':mref_ids['Recal_table_2'],'recal_table_0':mref_ids['Recal_table_0']}
            before_after_ids[before_after] = connection.add(package+'BQSR_'+before_after+'_grp', data)[0]
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id),
                'sh_script':sh_id,'bQSR_before_grp':before_after_ids['before'],'bQSR_after_grp':before_after_ids['after']}
        added_id = connection.add(package+'BQSR', data)[0]
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'bQSR':added_id})
def parse_addOrReplaceReadGroups(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start addOrReplaceGroups')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'AddOrReplaceReadGroups*.sh'), connection,package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'AddOrReplaceReadGroups', data)[0]
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'addOrReplaceReadGroups':added_id})
def parse_samToFilteredBam(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start samToFilteredBam')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'SamToFilteredBam*.sh'), connection,package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'SamToFilteredBam', data)[0]
        samToFilteredBam_data = connection.get(package+'SamToFilteredBam', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(samToFilteredBam_data) >0 and len(samToFilteredBam_data[0]['id']) > 0:
            added_id = samToFilteredBam_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'samToFilteredBam':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_sortBam(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print ('start sortBam')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'SortBam*.sh'), connection,package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'SortBam', data)[0]
        samToFilterdBam_data = connection.get(package+'SortBam', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(samToFilterdBam_data) >0 and len(samToFilterdBam_data[0]['id']) > 0:
            added_id = samToFilterdBam_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'sortBam':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_indelReallignmentKnown(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start indelReallignmentKnown')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'IndelReallignmentKnown*.sh'), connection,package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'IndelReallignmentKnow', data)   
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'indelReallignmentKnown':added_id})
def parse_gatkSplitNTrim(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start gatkSplitAndtrim')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'GATKSplitNTrim*.sh'), connection,package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'GATKSplitNTrim', data)[0]
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'gATKSplitNTrim':added_id})
def parse_cmMetrics(runinfo_folder,connection,package, pipeline):
    print ('start cmmMetrics')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder,'CollectMultipleMetrics*.sh'),connection, package):
        cmMetrics_base_file = re.search('O=(\S+)',sh_text).group(1)
        with open(cmMetrics_base_file+'.alignment_summary_metrics') as alignment_summary_metrics:
            metrics_class = alignment_summary_metrics.read().split('## METRICS CLASS')[1].split('\n')
            alignment_summary_metrics_ids = ''
            to_add = []
            for line in metrics_class[2:]:
                if len(line.strip()) == 0: continue
                s_l = line.split('\t')
                data = {'category':s_l[0],'total_reads':s_l[1],'pf_reads':s_l[2],'pct_pf_reads':s_l[3],'pf_noise_reads':s_l[4],
                        'pf_reads_aligned':s_l[5],'pct_pf_reads_aligned':s_l[6],'pf_aligned_bases':s_l[7],'pf_hq_aligned_reads':s_l[8],
                        'pf_h1_aligned_bases':s_l[9],'pf_hq_aligned_q20_bases':s_l[10],'pf_h1_median_mismatches':s_l[11],'pf_mismatch_rate':s_l[12],
                        'pf_hq_error_rate':s_l[13],'pf_indel_rate':s_l[14],'mean_read_length':s_l[15],'reads_aligned_in_pairs':s_l[16],
                        'pct_reads_aligned_in_pairs':s_l[17],'bad_cycles':s_l[18],'strand_balance':s_l[19],'pct_chimeras':s_l[20],
                        'pct_adapter':s_l[21],'sample':s_l[22],'library':s_l[23],'read_group':s_l[24]}
                to_add.append(data)
            alignment_summary_metrics = ','.join(add_multiple_rows(package+'Alignment_summary_metrics',to_add,connection,True))
        with open(cmMetrics_base_file+'.insert_size_metrics') as insert_size_metrics:
            split_file = insert_size_metrics.read().split('## METRICS CLASS')[1].split('## HISTOGRAM')[0]
            metrics_class = split_file[0].split('\n')
            insert_size_metrics_ids = ''
            to_add = []
            for line in metrics_class[2:]:
                if len(line.strip()) == 0: continue
                s_l = line.split('\t')
                data = {'median_insert_size':s_l[0],'median_absolute_deviation':s_l[1],'min_insert_size':s_l[2],'max_insert_size':s_l[3],
                        'mean_insert_size':s_l[4],'standard':s_l[5],'_deviation':s_l[6],'read_pairs':s_l[7],'pair_orientation':s_l[8],
                        'width_of_10_percent':s_l[9],'width_of_20_percent':s_l[10],'width_of_30_perc':s_l[11],'ent':s_l[12],
                        'width_of_40_percent':s_l[13],'width_of_50_percent':s_l[14],'width_of_60_percent':s_l[15],
                        'width_of_70_percent':s_l[16],'width_of_80_perc':s_l[17],'ent':s_l[18],'width_of_90_percent':s_l[19],
                        'width_of_99_percent':s_l[20],'sample':s_l[21],'library':s_l[22],'read_group':s_l[23]}
                to_add.append(data)
            insert_size_metrics_ids = ','.join(add_multiple_rows(package+'Insert_size_metrics_class',to_add,connection,True))
            histogram = split_file[1].split('\n')
            insert_size_histogram_ids = ''
            for line in histogram[2:]:
                if len(line.strip()) == 0: continue
                s_l = line.split('\t')
                data = {'insert_size':s_l[0],'all_reads_fr_count':s_l[1],'all_reads_rf_count':s_l[2]}
                insert_size_histogram_ids += connection.add(package+'Insert_size_metrics_histogram',data)+','
        with open(cmMetrics_base_file+'.quality_by_cycle_metrics') as alignment_summary_metrics:
            metrics_class = alignment_summary_metrics.read().split('## HISTOGRAM')[1].split('\n')
            quality_per_cycle_histogram_ids = ''
            to_add = []
            for line in metrics_class[2:]:
                if len(line.strip()) == 0: continue
                s_l = line.split('\t')
                data = {'cycle':s_l[0],'mean_quality':s_l[1]}
                to_add.append(data)
            quality_per_cycle_histogram_ids = ','.join(add_multiple_rows(package+'Quality_by_cycle_metrics',to_add,connection,True))                         
        with open(cmMetrics_base_file+'.quality_by_cycle_metrics') as alignment_summary_metrics:
            metrics_class = alignment_summary_metrics.read().split('## HISTOGRAM')[1].split('\n')
            quality_distribution_histogram_ids = ''
            to_add = []
            for line in metrics_class[2:]:
                if len(line.strip()) == 0: continue
                s_l = line.split('\t')
                data = {'cycle':s_l[0],'mean_quality':s_l[1]}
                to_add.append(data)
            quality_distribution_histogram_ids = ','.join(add_multiple_rows(package+'Quality_distribution_metrics',to_add,connection,True))  
        data = {'alignment_summary_metrics':alignment_summary_metrics_ids.rstrip(','),'insert_size_metrics_class':insert_size_metrics_ids.rstrip(','),'insert_size_metrics_histogram':insert_size_histogram_ids.rstrip(','),
                'pipeline':pipeline,'qual_by_cycle_metrics':quality_per_cycle_histogram_ids.rstrip(','),'qual_distribution_metrics':quality_distribution_histogram_ids.rstrip(','),
                'err_file':err_id,'out_file':out_id,'runtime':runtime,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'CMMetrics', data)[0]  
        if 'QC' in runinfo_folder:
            cMMetric_data = connection.get(package+'CMMetrics', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
            if len(cMMetric_data) >0 and len(cMMetric_data[0]['id']) > 0:
                added_id = cMMetric_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'cMMetrics':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_rMetrics(runinfo_folder,connection,package,pipeline):
    '''finished'''
    print ('start rMetrics')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder,'CollectRnaSeqMetrics*.sh'),connection, package):
        rMetrics = re.search('OUTPUT=(\S+.rna_metrics.log)', sh_text).group(1)
        # remove below line when on cluster
        with open(rMetrics) as rMetrics_file:
            split_file = rMetrics_file.read().split('## METRICS CLASS')[1].split('## HISTOGRAM')[0]
            metrics_class = split_file[0].split('\n')
            metrics_ids = ''
            to_add = []
            for line in metrics_class[2:]:
                if len(line.strip()) == 0: continue
                s_l = line.split('\t')
                data = {'pf_bases':s_l[0],'pf_aligned_bases':s_l[1],'ribosomal_bases':s_l[2],'coding_bases':s_l[3],
                        'utr_bases':s_l[4],'intronic_bases':s_l[5],'intergenic_bases':s_l[6],'ignored_reads':s_l[7],
                        'correct_strand_reads':s_l[8],'incorrect_strand_reads':s_l[9],'pct_ribosomal_bases':s_l[10],
                        'pct_coding_bases':s_l[11],'pct_utr_bases':s_l[12],'pct_intronic_bases':s_l[13],'pct_intergenic_bases':s_l[14],
                        'pct_mrna_bases':s_l[15],'pct_usable_bases':s_l[16],'pct_correct_strand_reads':s_l[17],
                        'median_cv_coverage':s_l[18],'median_5prime_bias':s_l[19],'median_3prime_bias':s_l[20],
                        'median_5prime_to_3prime_bias':s_l[21],'sample':s_l[22],'library':s_l[23],'read_group':s_l[24]}
                to_add.append(data)
            metrics_ids = ','.join(add_multiple_rows(package+'Rnaseq_metrics_class',to_add,connection,True))
            histogram = split_file[1].split('\n')
            metrics_histogram_ids = ''
            for line in histogram[2:]:
                if len(line.strip()) == 0: continue
                s_l = line.split('\t')
                data = {'all_reads_normalized_coverage':s_l[0],'normalized_position':s_l[1],'sample_exp_normalized_coverage':s_l[2],'sample_normalized_coverage':s_l[3]}
                metrics_histogram_ids += connection.add(package+'Rnaseq_metrics_histogram',data)+','
 
        data = {'rnaseq_metrics_class':metrics_ids.rstrip(','),'rnaseq_metrics_histogram':metrics_histogram_ids.rstrip(','),'err_file':err_id,'out_file':out_id,'runtime':runtime,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'CRMetrics', data)[0]
        if 'QC' in runinfo_folder:
            cRMetric_data = connection.get(package+'CRMetrics', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
            if len(cRMetric_data) >0 and len(cRMetric_data[0]['id']) > 0:
                added_id = cRMetric_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'cRMetrics':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_mergeBam(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start mergeBam')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'MergeBamFiles*.sh'),connection, package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,'internalId':internalId,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'MergeBamFiles', data)[0]
        rRMetric_data = connection.get(package+'rRMetric', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(rRMetric_data) >0 and len(rRMetric_data[0]['id']) > 0:
            added_id = rRMetric_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'rRMetric':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_flagstat(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start Flagstat')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'Flagstat*.sh'),connection, package):
        flagstat = re.search('if samtools flagstat \S+ \> (\S+.flagstat)',sh_text).group(1)
        # remove below line when on cluster
        #flagstat ='/Users/Niek/UMCG/test/data/ATACseq/project/flagStat/A_S1_L001.flagstat'
        with open(flagstat) as flagstat_file:
            flagstat_data = flagstat_file.read()
            groups = re.search('(\d+) \+ 0 in total.*?'+\
                               '(\d+) \+ 0 secondary.*?'+\
                               '(\d+) \+ 0 supplementary.*?'+\
                               '(\d+) \+ 0 duplicates.*?'+\
                               '(\d+) \+ 0 mapped \((\d+.\d+)%.*?'+\
                               '(\d+) \+ 0 paired in sequencing.*?'+\
                               '(\d+) \+ 0 read1.*?'+\
                               '(\d+) \+ 0 read2.*?'+\
                               '(\d+) \+ 0 properly paired \((\d+.\d+)%.*?'+\
                               '(\d+) \+ 0 with itself and mate mapped.*?'+\
                               '(\d+) \+ 0 singletons \((\d+.\d+)%.*?'+\
                               '(\d+) \+ 0 with mate mapped to a different chr.*?'+\
                               '(\d+) \+ 0 with mate mapped to a different chr \(mapQ\>=(\d+)\)',flagstat_data,re.DOTALL)
        data = {'mate_mapped_to_diff_chr':groups.group(14),'mate_mapped_to_diff_chr_mq5':groups.group(15),'n_read1':groups.group(7),'duplicates':groups.group(4),
                'n_read2':groups.group(8),'paired_in_sequencing':groups.group(6),'properly_paired':groups.group(9),'properly_paired_perc':groups.group(10),'secondary_mapped_reads':groups.group(2),'total_reads':groups.group(1),
                'singletons':groups.group(12),'singletons_perc':groups.group(13),'supplementary':groups.group(3),'total_mapped_reads':groups.group(5),'with_itself_and_mate_mapped':groups.group(11),
                'total_mapped_reads_perc':groups.group(5),'err_file':err_id,'out_file':out_id,'runtime':runtime,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'Flagstat', data)   
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'flagstat':added_id})
def parse_kallisto(runinfo_folder_quantification,connection,package):
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project, sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_quantification,'Kallisto*.sh'), connection, package):
        output_folder = re.search('kallistoDir="(\S+)"', sh_text).group(1)
        with open(output_folder+'run_info.json') as run_info_json:
            run_info_text = run_info_json.read()
            n_targets = re.search('"n_targets": (\d+)', run_info_text).group(1)
            n_bootstraps = re.search('“n_bootstraps”: (\d+)', run_info_text).group(1)
            index_version = re.search('"index_version": (\d+)', run_info_text).group(1)
        abundance = ''
        with open(output_folder+'abundance.tsv') as abundance_tsv:
            abundance_tsv.readline()
            to_add = []
            for line in abundance_tsv:
                spl_line = line.strip().split('\t')
                data = {'target_id':spl_line[0], 'length':spl_line[1],'eff_length':spl_line[2],
                        'est_counts':spl_line[3],'tpm':spl_line[4]}
                to_add.append(data)
            abundance = ','.join(add_multiple_rows(package+'Abundance', to_add,connection,True))
        abundance = abundance.rstrip(',')
        k_mer_length = re.search('k-mer length:\s+(\d+)', err_text).group(1)
        number_of_targets = re.search('number of targets:\s+(\S+)', err_text).group(1).replace(',','')
        number_of_kmers = re.search('number of k-mers:\s+(\S+)', err_text).group(1).replace(',','')
        n_of_equivalence_classes = re.search('number of equivalence classes:\s+(\S+)', err_text).group(1).replace(',','')
        running_mode = re.search('running in (.*?) mode', err_text).group(1)
        reads_processed = re.search('processed (\S+) reads', err_text).group(1).replace(',','')
        reads_pseudoaligned = re.search('(\S+) reads pseudoaligned', err_text).group(1).replace(',','')
        rounds_of_em = re.search('Expectation-Maximization algorithm ran for (\S+) rounds', err_text).group(1).replace(',','')
        data = {'k_mer_length':k_mer_length,'number_of_targets':number_of_targets,'number_of_kmers':number_of_kmers,'n_of_equivalence_classes':n_of_equivalence_classes,
                'running_mode':running_mode,'reads_processed':reads_processed,'reads_pseudoaligned':reads_pseudoaligned,'rounds_of_EM':rounds_of_em,
                'n_targets':n_targets,'n_bootstraps':n_bootstraps,'index_version':index_version,'abundance':abundance,
                'err_file':err_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id),
                'out_file':out_id,'runtime':runtime,'sh_script':sh_id}
        added_id = connection.add(package+'CombineBedFiles', data)   
        kallisto_data = connection.get(package+'Kallisto', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(kallisto_data) >0 and len(kallisto_data[0]['id']) > 0:
            added_id = kallisto_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'kallisto':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_combineBedFiles(runinfo_folder_QC,connection,package):
    '''finished'''
    print('start combineBedFiles')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_QC,'CombineBedFiles*.sh'), connection,package):
        combineBed_base = re.search('--out \S+(combinedFiles)',sh_text).group(1)
        # remove below line when on cluster
        #combineBed_base ='/Users/Niek/UMCG/test/data/ATACseq/project/combinedBED/combinedFiles'
        with open(combineBed_base+'.log') as combineBed_log_file:
            combineBed_log = combineBed_log_file.read()
            groups = re.search('(\d+) markers to be included.*?'+\
                               'Reading pedigree information from.*?'+\
                               '(\d+) individuals read from.*?'+\
                               '(\d+) individuals with nonmissing phenotypes.*?'+\
                               '(\d+) cases, (\d+) controls and (\d+) missing.*?'+\
                               '(\d+) males, (\d+) females, and (\d+) of unspecified sex.*?',combineBed_log,re.DOTALL)
            
        data = {'mssnip_list':'NotImplemented','n_cases':groups.group(4),'n_controls':groups.group(5),'n_females':groups.group(8),
                'n_males':groups.group(7),'n_markers_included':groups.group(1),'n_missing':groups.group(6),
                'n_nonmissing':groups.group(3),'n_nosex':groups.group(9),'n_of_mssnip':'NotImplemented','n_individual_read':groups.group(2),
                'err_file':err_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id),
                'out_file':out_id,'runtime':runtime,'sh_script':sh_id}
        added_id = connection.add(package+'CombineBedFiles', data)[0]
        combineBedFiles_data = connection.get(package+'CombineBedFiles', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(combineBedFiles_data) >0 and len(combineBedFiles_data[0]['id']) > 0:
            added_id = combineBedFiles_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'combineBedFiles':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def analyseCovariates(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start analyseCovariates')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id,err_id,out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'AnalyseCovariates*.sh'), connection,package):
        ac_intermediate = re.search('-csv (\S+.csv)',sh_text).group(1)
        with open(ac_intermediate) as ac_intermediate_file:
            ac_intermediate_file.readline()
            ac_intermediate_ids = ''
            to_add = []
            for line in ac_intermediate_file:
                s_l = line.split(',')
                data = {'covariate_value':s_l[0],'covariate_name':s_l[1],'event_type':s_l[2],'observations':s_l[3],'errors':s_l[4],
                        'empirical_quality':s_l[5],'average_reported_quality':s_l[6],'accuracy':s_l[7],'recalibration':s_l[8], 
                        'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
                to_add.append(data)
            ac_intermediate_ids = ','.join(add_multiple_rows(package+'CovariateAnalysisTable',to_add,connection,True))
                
            ac_intermediate_ids = ac_intermediate_ids.rstrip(',')
            
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,
                'sh_script':sh_id,'covariateAnalysisTable':ac_intermediate_ids, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'AnalyseCovariates', data)[0]   
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'analyseCovariates':added_id})
def parse_mergeGvcf(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print ('start mergeGvcf')
    for sh_text, err_text, mout_text, runtime, sample_name, internalId, project,sh_id,err_id,out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'MergeGvcf*.sh'),connection, package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'MergeGvcf', data)[0]
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'mergeGvcf':added_id})
def parse_genotypeHarmonizer(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start genotypeHarmonizer')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'GenotypeHarmonizer*.sh'),connection, package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name),'internalId':internalId,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'GenotypeHarmonizer', data)[0]
        genotypeHarmonizer_data = connection.get(package+'GenotypeHarmonizer', [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}])
        if len(genotypeHarmonizer_data) >0 and len(genotypeHarmonizer_data[0]['id']) > 0:
            added_id = genotypeHarmonizer_data[0]['id']+','+added_id
        connection.update_entity_rows(package+'Samples', data={'genotypeHarmonizer':added_id}, row_id = str(project)+'-'+str(sample_name)+'-'+str(analysis_id))
def parse_indelRealignmentKnown(runinfo_folder_genotypeCalling,connection,package):
    '''finished'''
    print('start indelRealignmentKnown')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project,sh_id, err_id, out_id, tool_ids in parse_rnaseq_tools(os.path.join(runinfo_folder_genotypeCalling,'IndelRealignmentKnown*.sh'),connection, package):
        data = {'err_file':err_id,'out_file':out_id,'runtime':runtime,
                'sh_script':sh_id, 'tools':tool_ids,'sample_id':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}
        added_id = connection.add(package+'IndelRealignmentKnown', data)[0] 
        connection.update_entity_rows(package+'Samples', query_list = [{'field':'id','operator':'EQUALS','value':str(project)+'-'+str(sample_name)+'-'+str(analysis_id)}], data = {'indelRealignmentKnown':added_id})

if __name__ == "__main__":
    connection = molgenis.Connect_Molgenis('http://localhost:8080',new_pass_file=False)
    runinfo_QC = configSectionMap("paths")['runinfo_folder_qc']
    runinfo_genotypeCalling = configSectionMap("paths")['runinfo_folder_genotypecalling']
    runinfo_folder_quantification = configSectionMap("paths")['runinfo_folder_quantification']
    package = configSectionMap("settings")['package']
    parse_mergeGvcf(runinfo_genotypeCalling,connection,package)
    parse_genotypeHarmonizer(runinfo_genotypeCalling,connection,package)
    parse_indelRealignmentKnown(runinfo_genotypeCalling,connection,package)
    
    #connection.delete_all_entity_rows(package+'_Samples')
    
