


import csv, os, argparse

parser = argparse.ArgumentParser(description='imput file with all project names which can me removed from prm storage')
parser.add_argument('--inn', help='file with projects which can be removed')
parser.add_argument('--out', help='sh script with rm file_to_be_removed.fq')
parser.add_argument('--prm', help='choose between prm02 or prm03')
args = parser.parse_args()

print('inn:', args.inn)
print ('out:', args.out)
print ('prm:', args.prm)



## makes list of all projects which can be removed
projects_to_remove_list = []

with open(args.inn) as inputfile:
    removeproject = csv.reader(inputfile)
    for i in removeproject:
        print (i[0])
        projects_to_remove_list.append(i[0]) #projects_to_remove_list.append(i[0]+'.csv')


# #takes information from sample sheet in project/run01/jobs/ folder
list_discarded = []
   
with open(args.out, 'w') as writefile:
    for project in projects_to_remove_list:
        #print(project)
        projectpath = '/groups/umcg-gaf/'+args.prm+'/projects/'+project+'/'
        dirsproject = os.listdir(projectpath)
        for files in dirsproject:
            #print (files)
            if os.path.isfile(str('/groups/umcg-gaf/'+args.prm+'/projects/'+project+'/'+files+'/results/'+project+'.csv')):
                with open('/groups/umcg-gaf/'+args.prm+'/projects/'+project+'/'+files+'/results/'+project+'.csv', 'r') as samplesheet:
                    print(samplesheet)
                    reader = csv.DictReader(samplesheet)
                    for row in reader:
                        sample_fq = ('\n'+'rm '+'/groups/umcg-gaf/'+args.prm+'/rawdata/ngs/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'_L'+row['lane']+'_'+row['barcode']+'*')
                        writefile.write(sample_fq)
                        discarded_fq = ('\n'+'rm '+'/groups/umcg-gaf/'+args.prm+'/rawdata/ngs/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'_L'+row['lane']+'_DISCARDED*')
                        demultiplex_log = ('\n'+'rm '+'/groups/umcg-gaf/'+args.prm+'/rawdata/ngs/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'_L'+row['lane']+'.demultiplex.log')
                        read_count_checks_for_pairs = ('\n'+'rm '+'/groups/umcg-gaf/'+args.prm+'/rawdata/ngs/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'/'+row['sequencingStartDate']+'_'+row['sequencer']+'_'+row['run']+'_'+row['flowcell']+'_L'+row['lane']+'.read_count_check_for_pairs.passed')
                        if discarded_fq in list_discarded:
                            next
                        else:
                            writefile.write(discarded_fq)
                            writefile.write(read_count_checks_for_pairs)
                            writefile.write(demultiplex_log)
                            list_discarded.append(discarded_fq)
                writefile.write('\n'+'mv -f '+'/groups/umcg-gaf/'+args.prm+'/projects/'+project+'/'+files+'/results/'+project+'.csv'+' '+'/groups/umcg-gaf/'+args.prm+'/projects/removedProjects/' )
                writefile.write('\n'+'rm -rf '+'/groups/umcg-gaf/'+args.prm+'/projects/'+project+'/')
            else:
                print ( '/groups/umcg-gaf/'+args.prm+'/projects/'+project+'/'+files+'/results/'+project+'.csv not found')
                continue





