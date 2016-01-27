'''
Created on Jul 31, 2015

@author: Niek
'''

import time
import argparse
import random
try:
    # Python 3
    import configparser
except ImportError:
    # Python 2.7
    import ConfigParser
from RNAseqParser import molgenis_wrapper
import sys
import requests

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

config = configparser.RawConfigParser()
config.read(r'PublicRNAseqParser/CONFIG')
def configSectionMap(section):
    configs = {}
    options = config.options(section)
    for option in options:
        configs[option] = config.get(section, option)
        if option == 'package' and not configs[option].endswith('_'):
            configs[option] += '_'
        if configs[option] == -1:
            print(("skip: %s" % option))
    return configs

parser = argparse.ArgumentParser(prog='RNAseq pipeline output parser',
                                 description='Command line interface for filling a Molgenis database with'+\
                                             'data from RNAseq analaysis tools used in an RNAseq analysis pipeline',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
# all the available options
parser.add_argument("--all", help="Parse all tools", action='store_true')
parser.add_argument("-a","--analyseCovariates", help="Parse analyseCovariates", action='store_true')
parser.add_argument("-b","--bqsr", help="Parse bqsr", action='store_true')
parser.add_argument("-c","--cMetrics_QC", help="Parse collectMultipleMetrics QC pipeline", action='store_true')
parser.add_argument("-d","--analyse_covariates", help="Parse AnalyseCovariates for pipeline", action='store_true')
parser.add_argument("-e","--ena", help="Parse ENA samples info", action='store_true')
parser.add_argument("-f","--mergeGvcf", help="Parse mergeGvcf", action='store_true')
parser.add_argument("-g","--gatkSplitNtrim", help="Parse gatkSplitNtrim", action='store_true')
parser.add_argument("-i","--indelRealignmentKnown", help="Parse indelRealignmentKnown", action='store_true')
parser.add_argument("-j","--genotypeHarmonizer", help="Parse genotypeHarmonizer", action='store_true')
parser.add_argument("-k","--markDuplicates", help="Parse markDuplicates", action='store_true')
parser.add_argument("-l","--mergeBam", help="Parse mergeBam", action='store_true')
parser.add_argument("-m","--sortBam", help="Parse sortBam", action='store_true')
parser.add_argument("-n","--flagstat", help="Parse flagstat", action='store_true')
parser.add_argument("-o","--addOrReplaceReadGroups", help="Parse addOrReplaceReadGroups", action='store_true')
parser.add_argument("-p","--combineBed", help="Parse combineBed", action='store_true')
parser.add_argument("-q","--fastqc", help="Parse fastqc", action='store_true')
parser.add_argument("-r","--rMetrics_QC", help="Parse RnaMetrics for QC pipeline", action='store_true')
parser.add_argument("-s","--samToFilteredBam", help="Parse samToFilteredBam", action='store_true')
parser.add_argument("-t","--hisat", help="Parse Hisat output", action='store_true')
parser.add_argument("-u","--unifiedGenotyper", help="Parse unifiedGenotyper", action='store_true')
parser.add_argument("-v","--variantEval", help="Parse variantEval", action='store_true')
parser.add_argument("-w","--haplotypeCaller", help="Parse haplotypeCaller", action='store_true')
parser.add_argument("-x","--cMetrics_genotypeCalling", help="Parse collectMultipleMetrics genotypeCalling pipeline", action='store_true')
parser.add_argument("-y","--verifyBamID", help="Parse verifyBamID", action='store_true')
parser.add_argument("-z","--rMetrics_genotypeCalling", help="Parse RnaMetrics for genotypeCalling pipeline", action='store_true')
parser.add_argument("-5", '--md5sum', help="parse md5sum", action='store_true')
parser.add_argument("-1", '--kallisto', help="parse Kallisto", action='store_true')
parser.add_argument("-2", '--gvcf', help="parse gvcf", action='store_true')
parser.add_argument("--QC", help="Add data from all QC steps",action='store_true')
parser.add_argument("--genotypeCalling", help="Add data from all genotypeCalling steps",action='store_true')
parser.add_argument("--delete_entity", help="Delete all rows of entity (happens before filling database)")
parser.add_argument("--delete_all", help="Delete all rows of all entites of package (happens before filling database)", action='store_true')
parser.add_argument("--analysis_id", help="Overwrite current analysis ID in CONFIG",default=configSectionMap("settings")['analysis_id'])
parser.add_argument("--runinfo_folder_qc", help="Overwrite runinfo_folder_qc in CONFIG",default=configSectionMap("paths")['runinfo_folder_qc'])
parser.add_argument("--runinfo_folder_genotypecalling", help="Overwrite runinfo_folder_genotypecalling in CONFIG",default=configSectionMap("paths")['runinfo_folder_genotypecalling'])
parser.add_argument("--runinfo_folder_quantification", help="Overwrite runinfo_folder_quantification in CONFIG",default=configSectionMap("paths")['runinfo_folder_quantification'])
parser.add_argument("--samplesheet", help="Overwrite samplesheet in CONFIG", default = configSectionMap("paths")['samplesheet'])
parser.add_argument("--analysis_description", help="Overwrite analysis_description in CONFIG",default = configSectionMap("settings")['analysis_description'])
parser.add_argument("--package", help="Overwrite package in CONFIG", default=configSectionMap('settings')['package'])
parser.add_argument("--server", help="Overwrite server in CONFIG", default=configSectionMap('settings')['server'])
parser.add_argument("--new_pass", help="Use a new password even if file with password is saved", action='store_true', default=configSectionMap('settings')['new_pass_file'])
parser.add_argument("--remove_pass_file", help="Remove the saved file with password for server", action='store_true', default=configSectionMap('settings')['remove_pass_file'])
parser.add_argument("--password_location", help="Location to save password to (save it to a folder only you have access to. If not possible, make sure remove_pass_file is set to True in Config file)", default=configSectionMap('settings')['password_location'])
parser.add_argument("--experiment_type", help="Change the experiment type in the Config file (rna_seq,atac_seq,pro_seq,unkown,hic-seq,dnase-seq,gro_seq)", default=configSectionMap("settings")['experiment_type'])
parser.add_argument("--ENA_path", help="Change the path to ENA info file in the Config file", default=configSectionMap("paths")['ena'])
parser.add_argument("--max_rows", help="Set the maximum amount of rows to be added at the same time", default=configSectionMap("settings")['max_rows'])
parser.add_argument("--project_folder", help="Set the project folder location", default=configSectionMap("paths")['project_folder'])
args = parser.parse_args()

# overwrite the values in the config file with any values given on the command line. 
# If none given, the overwritten value is same as current value in the config file
with open(r'PublicRNAseqParser/CONFIG','w') as configfile:
    config.set('paths','ena',args.ENA_path)
    config.set('settings','analysis_id',args.analysis_id)
    config.set('paths','runinfo_folder_qc',args.runinfo_folder_qc)
    config.set('paths','runinfo_folder_genotypecalling',args.runinfo_folder_genotypecalling)
    config.set('paths','samplesheet',args.samplesheet)
    config.set('settings','analysis_description',args.analysis_description)
    config.set('settings','package',args.package)
    config.set('settings','server',args.server)
    config.set('settings','new_pass_file',args.new_pass)
    config.set('settings','remove_pass_file',args.remove_pass_file)
    config.set('settings','password_location',args.password_location)
    config.set('settings','experiment_type',args.experiment_type)
    config.set('settings','max_rows', args.max_rows)
    config.set('settings','project_folder',args.project_folder)
    config.write(configfile)    
    
# This is imported here instead of at the top because if --max_rows is used 
# it won't be set until the second time you run it, as the input file will already be read
from PublicRNAseqParser import parse_output

print('Running parse_RNAseq_parser with configuration options:')
print((open('PublicRNAseqParser/CONFIG').read()))
  
# make a connection to the molgenis database. This connection will be passed to the other functions
with molgenis_wrapper.Connect_Molgenis(configSectionMap('settings')['server'],
                                remove_pass_file = configSectionMap('settings')['remove_pass_file'],
                                new_pass_file = configSectionMap('settings')['new_pass_file']) as connection:
    connection._add_datetime_default = True
    connection._added_by_default = True
    connection._updated_by_default = True
    print('connection established')
    try:
        rundir_QC = configSectionMap('paths')['runinfo_folder_qc']
        rundir_genotypeCalling = configSectionMap('paths')['runinfo_folder_genotypecalling']
        rundir_quantification = configSectionMap('paths')['runinfo_folder_quantification']
        package = configSectionMap('settings')['package']
    except KeyError:
        print('One of the runinfo paths missing in the CONFIG file')
        raise
    if args.all:
        # add ena later
        args.hisat, args.variantEval, args.verifyBamID, args.bqsr, args.addOrReplaceReadGroups, args.combineBed = (True,)*6
        args.samToFilteredBam, args.sortBam, args.indelRealignmentKnown, args.gatkSplitNtrim, args.rMetrics_QC, args.flagstat = (True,)*6
        args.cMetrics_QC, args.rMetrics_genotypeCalling, args.cMetrics_genotypeCalling, args.analyse_covariates, args.mergeBam = (True,)*5
        args.haplotypeCaller, args.unifiedGenotyper, args.genotypeHarmonizer, args.mergeGvcf, args.fastqc, args.markDuplicates = (True,)*6
        args.analyseCovariates, args.md5sum, args.kallisto, args.gvcf = (True,)*4
    if args.QC:
        args.hisat, args.variantEval, args.verifyBamID, args.combineBed = (True,)*4
        args.samToFilteredBam, args.sortBam, args.rMetrics_QC,args.cMetrics_QC = (True,)*4
        args.unifiedGenotyper, args.genotypeHarmonizer, args.fastqc,args.md5sum = (True,)*4
    if args.genotypeCalling:
        args.bqsr, args.addOrReplaceReadGroups,args.indelRealignmentKnown,args.gatkSplitNtrim, args.flagstat  = (True,)*5
        args.rMetrics_genotypeCalling, args.cMetrics_genotypeCalling, args.analyse_covariates, args.mergeBam = (True,)*4
        args.haplotypeCaller, args.mergeGvcf,args.markDuplicates,args.analyseCovariates,args.md5sum,args.gvcf = (True,)*3
    if args.delete_all:
        entities_list = []
        with open(configSectionMap('paths')['entities']) as entities:
            for entity in entities:
                entities_list.append(entity)
            index = 0
            print('Deleting all rows of all enities in entities.txt....')
            while len(entities_list) > 0:
                if index >= len(entities_list):
                    index = 0
                try:
                    entity = entities_list[index]
                    connection.delete_all_entity_rows(configSectionMap('settings')['package']+entity)
                    entities_list.remove(entity)
                except Exception as e:
                    message = e.response.json()['errors'][0]['message']
                    if 'foreign key constraint fails' in message or 'Cannot delete entity because there are other entities referencing it. Delete these first' in message:
                        pass
                    else:
                        raise
                index += 1
    elif args.delete_entity:
        connection.delete_all_entity_rows(configSectionMap('settings')['package']+args.delete_entity)
    # always make sure analysis exists
    try:
        connection.add(package+'Analysis_info', data=[{'id':configSectionMap('settings')['analysis_id'], 'analysis_description':configSectionMap('settings')['analysis_description']}])
    except requests.exceptions.HTTPError as e:
        if 'Duplicate value' in str(e):
            pass
    start_time = time.time()
    print('\n'+'~'*10+'STARTING'+'~'*10+'\n')
    parse_output.parse_samples( configSectionMap('paths')['samplesheet'],connection,package=package,experiment_type=configSectionMap('settings')['experiment_type'])
    if args.hisat:
        parse_output.parse_hisat(rundir_QC,connection,package=package)
    if args.ena:
        parse_output.parse_ena( configSectionMap('paths')['ena'],connection, package=package,project=configSectionMap('settings')['project'])
    if args.verifyBamID:
        parse_output.parse_verifyBamID(rundir_QC,connection, package)
    if args.samToFilteredBam:
        parse_output.parse_samToFilteredBam(rundir_QC,connection, package)
    if args.sortBam:
        parse_output.parse_sortBam(rundir_QC,connection, package)
    if args.rMetrics_QC:
        parse_output.parse_rMetrics(rundir_QC,connection, package,'QC')
    if args.cMetrics_QC:
        parse_output.parse_cmMetrics(rundir_QC,connection, package,'QC')
    if args.addOrReplaceReadGroups:
        parse_output.parse_addOrReplaceReadGroups(rundir_genotypeCalling,connection, package)
    if args.gatkSplitNtrim:
        parse_output.parse_gatkSplitNTrim(rundir_genotypeCalling,connection, package)
    if args.indelRealignmentKnown:
        parse_output.parse_indelRealignmentKnown(rundir_genotypeCalling,connection, package)
    if args.rMetrics_genotypeCalling:
        parse_output.parse_rMetrics(rundir_genotypeCalling,connection, package,'GenotypeCalling')
    if args.cMetrics_genotypeCalling:
        parse_output.parse_cmMetrics(rundir_genotypeCalling,connection, package,'GenotypeCalling')
    if args.analyse_covariates:
        parse_output.analyseCovariates(rundir_genotypeCalling, connection, package)
    if args.haplotypeCaller:
        parse_output.parse_variantCaller('HaplotypeCaller', rundir_genotypeCalling, connection, package)
    if args.unifiedGenotyper:
        parse_output.parse_variantCaller('UnifiedGenotyper', rundir_QC, connection, package)
    if args.gvcf:
        parse_output.parse_variantCaller('GenotypeGvcf', rundir_genotypeCalling, connection, package)
    if args.mergeGvcf:
        parse_output.parse_mergeGvcf(rundir_genotypeCalling, connection, package)
    if args.genotypeHarmonizer:
        parse_output.parse_genotypeHarmonizer(rundir_genotypeCalling, connection, package)
    if args.fastqc:
        parse_output.parse_fastqc(rundir_QC, connection, package)
    if args.markDuplicates:
        parse_output.parse_markDuplicates(rundir_genotypeCalling, connection, package)
    if args.flagstat:
        parse_output.parse_flagstat(rundir_genotypeCalling, connection, package)
    if args.combineBed:
        parse_output.parse_combineBedFiles(rundir_QC, connection, package)
    if args.md5sum:
        parse_output.parse_md5sums(connection, package)
    if args.variantEval:
        parse_output.parse_variantEval(rundir_QC,connection,package=package)
    if args.bqsr:
        parse_output.parse_bqsr(rundir_genotypeCalling,connection, package)
    if args.kallisto:
        parse_output.parse_kallisto(rundir_quantification,connection, package)
        
seconds = int(time.time() - start_time)
m, s = divmod(seconds, 60)
h, m = divmod(m, 60)
print('\n'+'~'*10+'FINISHED in: '+"%d:%02d:%02d" % (h, m, s)+'~'*10)
