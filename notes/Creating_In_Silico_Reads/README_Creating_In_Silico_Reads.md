# Steps to auto-validate our Variant Calling pipeline by spiking-in simulated reads
Rough idea is to simulate reads from phiX with known SNPs, spike-in those reads at the start of our pipeline, and if successful our pipeline should find and QC-report those SNPs.

## Download and compile simulation software
We have chosen to use 'wgsim' (`[https://github.com/lh3/wgsim]`) to simulate reads.

Download the software
```bash
	wget [https://github.com/lh3/wgsim/archive/master.zip]
```
Unzip it, and compile:
```bash
	gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
```

## Download phiX ref genome
On [http://www.ncbi.nlm.nih.gov/nuccore/9626372?report=fasta]{:target="_blank"} you'll find "Enterobacteria phage phiX174 sensu lato, complete genome". Save its sequence in 'phiX174_ref.fasta'


## Append phiX to human reference genome
cat /gcc/resources//b37/indices//human_g1k_v37.fasta /gcc/groups/gcc/home/mdijkstra/development/InSilicoData/data/humanPhiX/phiX174_ref.fasta > /gcc/groups/gcc/home/mdijkstra/development/InSilicoData/data/humanPhiX/human_g1k_v37_phiX.fasta

### Derive files
A dict file:
```bash
	module load picard-tools/1.102
	java -jar $PICARD_HOME/CreateSequenceDictionary.jar R= human_g1k_v37_phiX.fasta O= human_g1k_v37_phiX.dict
```

A .amb, .ann, .bwt, .pac, and .sa file:
```bash
	module load bwa/0.7.10
	bwa index human_g1k_v37_phiX.fasta
```

A .fai file:
```bash
	module load samtools/0.1.19
	samtools faidx human_g1k_v37_phiX.fasta
```

### Bed files and interval lists
Create directory `mkdir -p ~/development/InSilicoData/data/humanPhiX/intervals` (NB please update path). In that directory create BED file `None_b37_human_g1k_v37_phiX.bed` with following content:
```bash
1       1       249250621       chr1    1000    +
2       1       243199373       chr2    1000    +
3       1       198022430       chr3    1000    +
4       1       191154276       chr4    1000    +
5       1       180915260       chr5    1000    +
6       1       171115067       chr6    1000    +
7       1       159138663       chr7    1000    +
8       1       146364022       chr8    1000    +
9       1       141213431       chr9    1000    +
10      1       135534747       chr10   1000    +
11      1       135006516       chr11   1000    +
12      1       133851895       chr12   1000    +
13      1       115169878       chr13   1000    +
14      1       107349540       chr14   1000    +
15      1       102531392       chr15   1000    +
16      1       90354753        chr16   1000    +
17      1       81195210        chr17   1000    +
18      1       78077248        chr18   1000    +
19      1       59128983        chr19   1000    +
20      1       63025520        chr20   1000    +
21      1       48129895        chr21   1000    +
22      1       51304566        chr22   1000    +
X       1       155270560       chrX    1000    +
Y       1       59373566        chrY    1000    +
MT      1       16569   chrMT   1000    +
NC_001422.1     1       5386    phiX174 1000    +
```

Next symlink both baits and exons to that file:
```bash
	ln -s None_b37_human_g1k_v37_phiX.bed None_exons_b37_human_g1k_v37_phiX.bed
	ln -s None_b37_human_g1k_v37_phiX.bed None_baits_b37_human_g1k_v37_phiX.bed
```

Create interval lists for baits and exons:
```bash
	perl /gcc/tools/scripts/create_interval_listV4.pl -Ref ../human_g1k_v37_phiX.dict -Exons None_exons_b37_human_g1k_v37_phiX -baits None_baits_b37_human_g1k_v37_phiX
```

Split BED files:
```bash
	INPUTDIR=~/development/InSilicoData/data/humanPhiX/intervals
	NAME=None_exons_b37_human_g1k_v37_phiX
	OUTPUTDIR=~/development/InSilicoData/data/humanPhiX/intervals

	if [ -f ${OUTPUTDIR}/${NAME}.chr1.bed ]
	        then
	        echo "Is this bed file already splitted before? If so, please remove the old ones or do not run this script ;)" 
	else
	        echo "splitting bed file"
	        awk '{print $0 >> "'${OUTPUTDIR}'/'${NAME}'.chr" $1".bed" }' ${INPUTDIR}/${NAME}.bed 
	        awk '{
	        if ($1 != "X"){
	                print $0 >> "'${OUTPUTDIR}'/'${NAME}'.withoutChrX.bed"
	        }' ${INPUTDIR}/${NAME}.bed
	fi
```

## Simulate the reads
```bash
	../wgsim/wgsim -N10000 -1 100 -2 100 phiX174_ref.fasta out1.fq out2.fq
```

Which gives output like:
```bash
	[wgsim] seed = 1430742799
	[wgsim_core] calculating the total length of the reference sequence...
	[wgsim_core] 1 sequences, total length: 5386
	gi|9626372|ref|NC_001422.1|	1574	G	R	+
	gi|9626372|ref|NC_001422.1|	1652	G	-	-
	gi|9626372|ref|NC_001422.1|	1915	T	C	-
```	

Please see [http://www.bioinformatics.org/sms/iupac.html] for substitution rules (e.g. R = A or G).

### Align reads and quickly check results
```bash
	bwa mem \
	-M \
	-R "@RG\tID:FakePhiX\tSM:FakePhiXSample\tPL:SimulatedByWgsim" \
	-t 8 \
	./phiX174_ref.fasta \
	out1.fq \
	out2.fq \
	> fakePhiX.sam
```
Check by `../wgsim/wgsim_eval.pl unique fakePhiX.sam | ../wgsim/wgsim_eval.pl alneval -g 20`
Additionally visually inspect by `samtools tview fakePhiX.sorted.bam phiX174_ref.fasta`.

### Standardize and prepare simulation files for analysis
Quality of all nucleotides simulated with wgsim have quality 2. You may 'improve' the quality, standardize file names, create md5 checksums and save in 'rawdata' directory by:

```bash
	awk ' /222222222/ { gsub("2", "I"); print $0; next } { print } ' out1.fq > 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_1.fq
	awk ' /222222222/ { gsub("2", "I"); print $0; next } { print } ' out2.fq > 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_2.fq

	gzip -c 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_1.fq > 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_1.fq.gz
	gzip -c 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_2.fq > 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_2.fq.gz

	md5sum 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_1.fq.gz > 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_1.fq.gz.md5
	md5sum 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_2.fq.gz > 150504_WGSIM_9999_ZZZZZZZZXX_L9_ZZZZZZ_2.fq.gz.md5
	
	sudo -u gaf mkdir /gcc/groups/gaf/prm02/rawdata/ngs/150504_WGSIM_9999_ZZZZZZZZXX
	sudo -u gaf cp *.gz* /gcc/groups/gaf/prm02/rawdata/ngs/150504_WGSIM_9999_ZZZZZZZZXX/
```

## Test our analysis pipeline
### Prepare for 1th stage
Create and go to directory
```bash
	mkdir /gcc/groups/gaf/tmp01/generatedscripts/InSilicoData
	cd /gcc/groups/gaf/tmp01/generatedscripts/InSilicoData
```

Save following sample sheet `InSilicoData.csv` there:
```bash
	internalSampleID,externalSampleID,project,sequencer,contact,validationLog,labStatusPhase,labStatusComments,lastUpDate,sequencingStartDate,run,flowcell,lane,barcodeMenu,seqType,prepKit,capturingKit,arrayFile,arrayID,GAF_QC_Name,GAF_QC_Date,GAF_QC_Status,GCC_Analysis,GCC_QC_Name,GCC_QC_Date,GCC_QC_Status,TargetDateShipment,DataShippedDate,DataShippedTo,DataShippedBy,Comments,barcode,barcodeType
	1,InSilicoSampleID0,InSilicoData,WGSIM,mdijkstra,,,,,150504,9999,ZZZZZZZZXX,9,NO_BARCODE,PE,None,None,,,,,,Yes,,,,,,,,,ZZZZZZ,SIM
	2,InSilicoSampleID1,InSilicoData,WGSIM,mdijkstra,,,,,150504,9999,ZZZZZZZZXX,9,NO_BARCODE,PE,None,None,,,,,,Yes,,,,,,,,,ZZZZZZ,SIM
```

Make copy of paramters file `cp /gcc/tools/NGS_DNA-2.1.0/parameters.csv .` and update references to index (NB please update exact directory):
```bash
	indexFileID,human_g1k_v37_phiX
	indexFileDictionary,/gcc/groups/gcc/home/mdijkstra/development/InSilicoData/data/humanPhiX/human_g1k_v37_phiX.dict
	indexFileIDtest,${indexFileID}
	indexFile,/gcc/groups/gcc/home/mdijkstra/development/InSilicoData/data/humanPhiX/human_g1k_v37_phiX.fasta
	intervalListDir,/gcc/groups/gcc/home/mdijkstra/development/InSilicoData/data/humanPhiX/intervals/
```

Create following script `generate_run_dna2.1.0.sh` to generate our pipeline:
```bash
	#!/bin/bash

	PROJECT=InSilicoData
	TMPDIR=tmp01
	RUNID=test_pbs02

	NGS_DNA_HOME=/gcc/tools/NGS_DNA-2.1.0/
	module load molgenis-compute/v5_20150211
	module list

	if [ -f .compute.properties ];
	then
	     rm .compute.properties
	fi

	if [ -f /gcc/groups/gaf/tmp01/generatedscripts/${PROJECT}/out.csv  ];
	then
	        rm -rf /gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/out.csv
	fi

	sh $MC_HOME/convert.sh /gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/parameters.csv \
	/gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/out.csv

	sh $MC_HOME/molgenis_compute.sh \
	-p /gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/out.csv \
	-p $NGS_DNA_HOME/chrParameters.csv \
	-p /gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/${PROJECT}.csv \
	-w $NGS_DNA_HOME/create_in-house_ngs_projects_workflow.csv \
	-rundir /gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/scripts \
	--runid ${RUNID} \
	-o "workflowpath=$NGS_DNA_HOME/workflow.csv;\
	outputdir=scripts/jobs;mainParameters=/gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/out.csv;\
	chrParameters=$NGS_DNA_HOME/chrParameters.csv;\
	worksheet=/gcc/groups/gaf/${TMPDIR}/generatedscripts/${PROJECT}/${PROJECT}.csv" \
	-weave \
	--generate
```
### Running our pipeline
In `/gcc/groups/gaf/tmp01/generatedscripts/InSilicoData/` run `sh generate_run_dna2.1.0.sh` to launch 1th-stage of our pipeline rocket, then `cd /gcc/groups/gaf/tmp01/generatedscripts/InSilicoData/scripts` and run `submit.sh` to launch 2th-stage of rocket. Please go to directory the analysis pipeline should be generated by now: `cd /gcc/groups/gaf/tmp01/projects/InSilicoData/test_pbs02/jobs/`, and shorten predicted runtimes so those jobs will be scheduled quickly!
```bash
	perl -pi -e 's|#PBS -l walltime=..:..:..|#PBS -l walltime=00:10:00|g' *.sh
```
Then start pipeline by `sh submit.sh` and monitor progress with `qstat -umdijkstra` or `ls -ltr *.finished`.
## Validate output
Step `s13_VariantCalling_23.sh` is the script that should find the SNPs that are spiked in the simulated phiX reads. According to that script, you should see the SNPs by
```bash
	tail -4 /gcc//groups/gaf/tmp01/tmp//InSilicoData/test_pbs02//InSilicoData.chrNC_001422.1.variant.calls.vcf
```
Here we spiked in only three variants. Please update the length of the tail dependent on the total number of variants.




























