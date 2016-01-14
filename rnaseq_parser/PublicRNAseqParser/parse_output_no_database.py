'''
Created on Jul 4, 2015

@author: Niek

TODO: Make sure all non-rnaseqtool tables also get duplicate value prevention
'''
import warnings
import os
import re
import glob
import datetime
import zipfile

       
def parse_rnaseq_tools(sh_file_path):
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

        yield sh_text, err_text, out_text, runtime, sample_name, internalId, project
def parse_fastqc(runinfo_folder_QC):
    '''not working yet, row too long'''
    print('Start fastqc')
    for sh_text, err_text, out_text, runtime, sample_name, internalId, project in parse_rnaseq_tools(os.path.join(runinfo_folder_QC,'FastQC*.sh')):
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
                    graphs[description] = img_data
                                                              
                if 'fastqc_data.txt' in name:
                    fastqc_data = archive.read(name)
        fastqc_groups = re.search('Filename\s+(\S+)\n'\
                                 +'File type\s+(.+)\n'\
                                 +'Encoding\s+(.+)\n'\
                                 +'Total Sequences\s+(\d+)\n'\
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
        
        fastqc_data = {'pbsq':[],'ptsq':[],'psqs':[],'pbsc':[],'psgc':[],'pbnc':[],'sld':[],'sdl':[],'os':[],'ac':[],'kc':[]}
        to_add = []
        for line in fastqc_groups.group(8).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'base_letter':s_l[0], 'mean':s_l[1], 'median':s_l[2], 'lower_quartile':s_l[3], 'upper_quartile':s_l[4],'percentile_10th':s_l[5], 'percentile_90th':s_l[6]})
        fastqc_data['pbsq'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(9).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'tile':s_l[0], 'base_letter':s_l[1], 'mean':s_l[2]})
        fastqc_data['ptsq'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(10).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'quality':s_l[0], 'count':s_l[1]})
        fastqc_data['psqs'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(11).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'base_letter':s_l[0], 'a':s_l[1], 'g':s_l[2], 'c':s_l[3], 't':s_l[4]})
        fastqc_data['pbsc'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(12).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'gc_content':s_l[0], 'count':s_l[1]})
        fastqc_data['psgc'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(13).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'base_letter':s_l[0], 'n_count':s_l[1]})
        fastqc_data['pbnc'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(14).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'length':s_l[0], 'count':s_l[1]})
        fastqc_data['sld'].append(to_add)
        total_deduplicated_perc = fastqc_groups.group(15).split('\n')[0].split('\t')[1]            
        to_add = []
        for line in fastqc_groups.group(15).split('\n')[2:]:
            s_l = line.split('\t')
            to_add.append({'duplication_level':s_l[0],'perc_of_deduplicated':s_l[1],'perc_of_total':s_l[2]})
        fastqc_data['sdl'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(16).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'sequence':s_l[0],'count':s_l[1],'percentage':s_l[2],'possible_source':s_l[3]})
        fastqc_data['os'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(17).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'position':s_l[0], 'illumina_universal_adapter':s_l[1], 'illumina_small_rna_adapter':s_l[2], 'nextera_transposase_seq':s_l[3], 'solid_small_rna_adapter':s_l[4]})
        fastqc_data['ac'].append(to_add)
        to_add = []
        for line in fastqc_groups.group(18).split('\n')[1:]:
            s_l = line.split('\t')
            to_add.append({'seq':s_l[0], 'count':s_l[1], 'p_value':s_l[2], 'obs_exp_max':s_l[3], 'max_obs_exp_position':s_l[4]})
        fastqc_data['kc'].append(to_add)

        data = {'adapter_content_graph':graphs['adapter_content'],'kmer_content_graph':graphs['kmer_profiles'],'per_base_n_content_graph':graphs['per_base_n_content'],
                'per_base_seq_content_graph':graphs['per_base_sequence_content'],'per_base_seq_qual_graph':graphs['per_base_quality'],'per_seq_GC_content_graph':graphs['per_sequence_gc_content'],
                'per_seq_qual_scores_graph':graphs['per_sequence_quality'],'per_tile_seq_qual_graph':graphs['per_tile_quality'],'seq_duplication_levels_graph':graphs['duplication_levels'],
                'seq_length_distribution_graph':graphs['sequence_length_distribution'],
                'total_deduplicated_perc':total_deduplicated_perc,'adapter_content':fastqc_data['ac'].rstrip(','),'adapter_content_check':image_data['Adapter Content'][1],'basic_statistics_check':image_data['Basic Statistics'][1],'gc_perc':fastqc_groups.group(7),
                'kmer_content':fastqc_data['kc'].rstrip(','),'kmer_content_check':image_data['Kmer Content'][1],'overrepresented_seqs':fastqc_data['os'].rstrip(','),
                'overrepresented_seqs_check':image_data['Overrepresented sequences'][1],'per_base_N_content':fastqc_data['pbnc'].rstrip(','),'per_base_n_content_check':image_data['Per base N content'][1],'per_base_seq_content':fastqc_data['pbsc'].rstrip(','),
                'per_base_seq_content_check':image_data['Per base sequence content'][1],'per_base_seq_qual':fastqc_data['pbsq'].rstrip(','),
                'per_base_seq_qual_check':image_data['Per base sequence quality'][1],'per_seq_GC_content':fastqc_data['psgc'].rstrip(','),'per_seq_gc_content_check':image_data['Per sequence GC content'][1],
                'per_seq_qual_scores':fastqc_data['psqs'].rstrip(','),'per_seq_qual_scores_check':image_data['Per sequence quality scores'][1],'per_tile_seq_qual':fastqc_data['ptsq'].rstrip(','),
                'per_tile_seq_qual_check':image_data['Per tile sequence quality'][1],'seq_duplication_levels':fastqc_data['sdl'].rstrip(','),
                'seq_duplication_levels_check':image_data['Sequence Duplication Levels'][1],'seq_length_min':fastqc_groups.group(6).split('-')[0],'seq_length_distribution':fastqc_data['sld'].rstrip(','),
                'seq_length_distribution_check':image_data['Sequence Length Distribution'][1],'seqs_flagged_as_poor':fastqc_groups.group(5),'total_seqs':fastqc_groups.group(4),'seq_length_max':fastqc_groups.group(6).split('-')[1],'file_name':fastqc_groups.group(2),'file_type':fastqc_groups.group(2),'encoding':fastqc_groups.group(3),
                'internalId':internalId,'internalId_sampleid':internalId+'_'+str(project)+'-'+str(sample_name)}
        

if __name__ == "__main__":
    parse_fastqc('')
    
