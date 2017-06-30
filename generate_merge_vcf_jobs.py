
import math
import glob  
import argparse

parser = argparse.ArgumentParser(description='Generate slurm jobs for merging vcfs.')
parser.add_argument(['-o','--outdir'], required=True, action='store_const',help="Output directory")
parser.add_argument(['-l','--list_of_vcfs'], required=True, action='store_const',help="File with list of vcfs to merge")
parser.add_argument(['-j','--jobs_dir'], required=True, action='store_const',help="Dir to write jobs to")
parser.add_argument(['-r','--ref_genome'], required=True, action='store_const',help="Reference genome")

args = parser.parse_args()


outdir=args.outdir
list_of_vcfs = args.list_of_vcfs
jobs_dir = args.jobdir
ref_genome = args.ref_genome

template = """#!/bin/bash
#SBATCH --job-name=MergeGvcfs_batchREPLACEBATCH_chrREPLACECHROMOSOME
#SBATCH --output=MergeGvcfs_batchREPLACEBATCH_chrREPLACECHROMOSOME.out
#SBATCH --error=MergeGvcfs_batchREPLACEBATCH_chrREPLACECHROMOSOME.err
#SBATCH --qos=leftover
#SBATCH --time=1-23:59:59
#SBATCH --cpus-per-task 3
#SBATCH --mem 18gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


ENVIRONMENT_DIR="."
set -e
set -u

echo "## "$(date)" Start $0"

#Load gatk module
module load GATK/3.4-0-Java-1.7.0_80
module list

mkdir -p REPLACEOUTDIR


if java -Xmx16g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${TMPDIR} \
    -jar ${EBROOTGATK}/GenomeAnalysisTK.jar \
    -T CombineGVCFs \
    -R REPLACEREFGENOME \
    -o REPLACEOUTDIR/REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz \
    -L /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.chrREPLACECHROMOSOME.interval_list \
    REPLACEINPUT

then
 echo "returncode: $?"; 

if [ ! -f REPLACEOUTDIR/REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz ]; then
    echo "REPLACEOUTDIR/REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz"
    exit 1
fi
cd REPLACEOUTDIR/
md5sum REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz > REPLACEPROJECT.batchREPLACEBATCH_chrREPLACECHROMOSOME.g.vcf.gz.md5
 cd -
 echo "succes moving files";
else
 echo "returncode: $?";
 echo "fail";
fi


touch REPLACEFINISHED

echo "## "$(date)" ##  $0 Done "


"""

with open(list_of_vcfs) as input:
    samples = 0
    lines = input.read().split('\n')
for line in lines:
    if len(line.strip()) == 0:
        continue
    samples += 1
batches = math.ceil(samples/200)
print('writing to',batches,'batches')
batch = 0
input = {'1':{},'2':{},'3':{},'4':{},'5':{},'6':{},'7':{},
         '8':{},'9':{},'10':{},'11':{},'12':{},'13':{},'14':{},
         '15':{},'16':{},'17':{},'18':{},'19':{},'20':{},
         '21':{},'22':{},'23':{},'24':{},'25':{}}

lines = sorted(lines)
for line in lines:
    if len(line.strip()) == 0:
        continue
    for chr in range(1, 26, 1):
        if str(batch) in input[str(chr)]:
            input[str(chr)][str(batch)] += ' --variant '+line.strip().replace('chr15','chr'+str(chr))
        else:
            input[str(chr)][str(batch)] = ' --variant '+line.strip().replace('chr15','chr'+str(chr))            
    batch += 1
    if batch >= batches:
        batch = 0

for chr in input:
    for batch in input[chr]:
        outfile = jobs_dir+'MergeGvcfs_batch'+batch+'_chr'+chr+'.sh'
        with open(outfile,'w') as out:
            new_template = template.replace('REPLACEBATCH',batch)
            new_template = new_template.replace('REPLACECHROMOSOME',chr)
            new_template = new_template.replace('REPLACEOUTDIR',outdir)
            new_template = new_template.replace('REPLACEINPUT', input[chr][batch])
            new_template = new_template.replace('REPLACEPROJECT','BIOS')
            new_template = new_template.replace('REPLACEREFGENOME',ref_genome)
            new_template = new_template.replace('REPLACEFINISHED','MergeGvcfs_batch'+batch+'_chr'+chr+'.sh.finished')
            out.write(new_template)
        print(outfile)

