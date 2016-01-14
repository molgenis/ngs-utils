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

       

def parse_fastqc(fastqc_folder):
    '''not working yet, row too long'''
    print('Start fastqc')
    for output in glob.glob(fastqc_folder+'/*zip'):
        print(output)
        fastqc_output_base = output.rstrip('.zip')
        try:
            with open(fastqc_output_base+'.html') as fastqc_html_file:
                fastqc_html = fastqc_html_file.read()
                images = re.findall('img src="data:image/png;base64,(\S+)" alt="(\S+)"/>\<a href="\S+"\>(.*?)\<',fastqc_html)
                image_data = {}
                for image in images:
                    image_data[image[2]] = [image[0],image[1].strip('[').strip(']').lower()]
        except:
            print('error')
        with zipfile.ZipFile(fastqc_output_base+'.zip') as archive:
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
                                 +'Sequence length\s+(\S+?)\n'\
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
                'total_deduplicated_perc':total_deduplicated_perc,'adapter_content':fastqc_data['ac'],'adapter_content_check':image_data['Adapter Content'][1],'basic_statistics_check':image_data['Basic Statistics'][1],'gc_perc':fastqc_groups.group(7),
                'kmer_content':fastqc_data['kc'],'kmer_content_check':image_data['Kmer Content'][1],'overrepresented_seqs':fastqc_data['os'],
                'overrepresented_seqs_check':image_data['Overrepresented sequences'][1],'per_base_N_content':fastqc_data['pbnc'],'per_base_n_content_check':image_data['Per base N content'][1],'per_base_seq_content':fastqc_data['pbsc'],
                'per_base_seq_content_check':image_data['Per base sequence content'][1],'per_base_seq_qual':fastqc_data['pbsq'],
                'per_base_seq_qual_check':image_data['Per base sequence quality'][1],'per_seq_GC_content':fastqc_data['psgc'],'per_seq_gc_content_check':image_data['Per sequence GC content'][1],
                'per_seq_qual_scores':fastqc_data['psqs'],'per_seq_qual_scores_check':image_data['Per sequence quality scores'][1],'per_tile_seq_qual':fastqc_data['ptsq'],
                'per_tile_seq_qual_check':image_data['Per tile sequence quality'][1],'seq_duplication_levels':fastqc_data['sdl'],
                'seq_duplication_levels_check':image_data['Sequence Duplication Levels'][1],'seq_length_min':fastqc_groups.group(6),'seq_length_distribution':fastqc_data['sld'],
                'seq_length_distribution_check':image_data['Sequence Length Distribution'][1],'seqs_flagged_as_poor':fastqc_groups.group(5),'total_seqs':fastqc_groups.group(4),'seq_length_max':fastqc_groups.group(6),'file_name':fastqc_groups.group(2),'file_type':fastqc_groups.group(2),'encoding':fastqc_groups.group(3)}
        print(data)

parse_fastqc('fastqc/')
    
