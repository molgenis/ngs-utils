import os
import sys
from optparse import OptionParser
import glob

beagleImputed = '/target/gpfs2/lifelines_rp/releases/LL3/BeagleImputedPedAndMap'
beagleDosage = '/target/gpfs2/lifelines_rp/releases/LL3/BeagleImputedDosage'

success = True
rep = ''
parser = OptionParser()
parser.add_option('--studyId', help='Required to set study id e.g.:  "<study>"')
parser.add_option('--studyDir', help='Required to set study dir e.g.: "/target/gpfs2/lifelines_rp/releases/<release>/lifelines_<study>"')
parser.add_option('--mappingFile', help='Required to set mapping file e.g. "/target/gpfs2/lifelines_rp/releases/<release>/<study>.txt"')
(options, args) = parser.parse_args()
if options.studyId == None or options.studyDir == None or options.mappingFile == None:
    parser.print_help()
    sys.exit()

if not os.path.exists(options.studyDir):
    print options.studyDir + ' not exists'
    sys.exit()
if not os.path.exists(options.mappingFile):
    print options.mappingFile + ' not exists'
    sys.exit()
if not os.path.exists(options.studyDir + '/test'):
    os.makedirs(options.studyDir + '/test')
    
if not os.path.exists(options.studyDir + '/test/original.ids'):
    cmd = 'awk \'{print $2}\' ' + options.mappingFile + ' | tr -d \'\015\' | sort > ' + options.studyDir + '/test/original.ids'
    os.system(cmd)
    
for item in os.listdir(options.studyDir):
    if os.path.splitext(item)[1] == '.dose':
        doseFilePath = options.studyDir + '/' + item
        doseIdsFilePath = options.studyDir + '/test/' + item + '.ids'
        doseSnpFilePath = options.studyDir + '/test/' + item + '.snp'
        if not os.path.exists(doseIdsFilePath):
            cmd = 'head -n 1 ' + doseFilePath + ' | sed \'s/SNP\\tA1\\tA2\\t1\\t//g\' | sed \'s/\\t1\\t/\\n/g\' | tr -d \'\015\' > ' + doseIdsFilePath
            print cmd
            os.system(cmd)             
        if not os.path.exists(doseSnpFilePath):
            cmd = 'awk \'{print $1}\' ' + doseFilePath + ' | sed 1d > ' + doseSnpFilePath
            print cmd
            os.system(cmd)        
    if os.path.splitext(item)[1] == '.ped':
        pedFilePath = options.studyDir + '/' + item
        pedIdsFilePath = options.studyDir + '/test/' + item + '.ids'
        if not os.path.exists(pedIdsFilePath):
            cmd = 'awk \'{print $2}\' ' + pedFilePath + ' > ' + pedIdsFilePath
            print cmd
            os.system(cmd)
    if os.path.splitext(item)[1] == '.map':
        mapFilePath = options.studyDir + '/' + item
        mapSnpFilePath = options.studyDir + '/test/' + item + '.snp'
        if not os.path.exists(mapSnpFilePath):
            cmd = 'awk \'{print $2}\' ' + mapFilePath + ' > ' + mapSnpFilePath
            print cmd
            os.system(cmd)              

for item in os.listdir(beagleDosage):
    if os.path.splitext(item)[1] == '.dose':
        biDoseFilePath = beagleDosage + '/' + item
        biDoseSnpFilePath = beagleDosage + '/' + item + '.snp'
        if not os.path.exists(biDoseSnpFilePath):
            cmd = 'awk \'{print $1}\' ' + biDoseFilePath + ' | sed 1d > ' + biDoseSnpFilePath
            print cmd
            os.system(cmd)

for item in os.listdir(beagleImputed):
    if os.path.splitext(item)[1] == '.map':
        biMapFilePath = beagleImputed + '/' + item
        biMapSnpFilePath = beagleImputed + '/' + item + '.snp'
        if not os.path.exists(biMapSnpFilePath):
            cmd = 'awk \'{print $2}\' ' + biMapFilePath + ' > ' + biMapSnpFilePath
            print cmd
            os.system(cmd)

def testZeroDiffLines(cmd):
    global success
    c = int(os.popen(cmd + ' | wc -l').read())
    if c != 0:
        success = False
        report('FAIL ' + str(c) + ' DIFFERENCE(S). COMMAND IS:')
        report(cmd)

def fail(msg):
    global success 
    success = False
    report(msg)
    testResult()

def testResult():
    global success
    global rep
    if success:
        report('\n-------\nTEST OK\n-------\n')
    else:
        report('\n-----------\nTEST FAILED\n-----------\n')
    fw = open(options.studyDir + '/test/' + options.studyId + '_TESTREPORT.TXT','w')
    fw.write(rep)
    sys.exit()

def report(msg):
    global rep
    rep += '\n' + msg
    print msg

report('\n-----------------------\nSTARTING REPORT OUTPUT\n-----------------------\n')
report('\nInput parameters: ' + str(options))
for n in range(1, 23):       
    doseSnpSrcFilePath = beagleDosage + '/ImputedGenotypeDosageFormatPLINK-Chr' + str(n) + '.dose.snp'
    mapSnpSrcFilePath = beagleImputed + '/output.' + str(n) + '.map.snp'
    doseSnpSubFilePath = options.studyDir + '/test/' + options.studyId + '_chr' + str(n) + '.dose.snp'
    mapSnpSubFilePath = options.studyDir + '/test/' + options.studyId + '_chr' + str(n) + '.map.snp'
    doseIdsFilePath = options.studyDir + '/test/' + options.studyId + '_chr' + str(n) + '.dose.ids'
    pedIdsFilePath = options.studyDir + '/test/' + options.studyId + '_chr' + str(n) + '.ped.ids'
    if not os.path.exists(doseSnpSrcFilePath):
        fail('FAILED MISSING: ' + doseSnpSrcFilePath)
    if not os.path.exists(mapSnpSrcFilePath):
        fail('FAILED MISSING: ' + mapSnpSrcFilePath)
    if not os.path.exists(doseSnpSubFilePath):
        fail('FAILED MISSING: ' + doseSnpSubFilePath)
    if not os.path.exists(mapSnpSubFilePath):
        fail('FAILED MISSING: ' + mapSnpSubFilePath)
    if not os.path.exists(doseIdsFilePath):
        fail('FAILED MISSING: ' + doseIdsFilePath)
    if not os.path.exists(pedIdsFilePath):
        fail('FAILED MISSING: ' + pedIdsFilePath)
    cmd = 'diff --side-by-side --suppress-common-lines ' + doseSnpSrcFilePath + ' ' + mapSnpSrcFilePath
    testZeroDiffLines(cmd)
    cmd = 'diff --side-by-side --suppress-common-lines ' + mapSnpSrcFilePath + ' ' + doseSnpSubFilePath
    testZeroDiffLines(cmd) 
    cmd = 'diff --side-by-side --suppress-common-lines ' + doseSnpSubFilePath + ' ' + mapSnpSubFilePath
    testZeroDiffLines(cmd)
    cmd = 'sort ' + doseIdsFilePath + ' | diff - --side-by-side --suppress-common-lines ' + options.studyDir + '/test/original.ids'
    testZeroDiffLines(cmd)
    cmd = 'sort ' + pedIdsFilePath + ' | diff - --side-by-side --suppress-common-lines ' + options.studyDir + '/test/original.ids'
    testZeroDiffLines(cmd)

report('\nNumber of participants: ' + os.popen('cat ' + options.studyDir + '/test/original.ids | wc -w').read())
testResult()
