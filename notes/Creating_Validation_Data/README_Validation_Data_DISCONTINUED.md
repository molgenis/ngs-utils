# Steps to get data for NGS pipeline validation
=== This work is currently discontinued. ===
=== Please see README_Validate_Delly_With_VarSim.md for a continuation Validating the INDEL calling step in our pipeline. ===

Currently, we want to validate our pipeline, c.q. its INDEL calling step, on data that originates from three different sources:
1. One, WGS set, which can be downloaded from http://www.ncbi.nlm.nih.gov/sra/?term=PRJNA281509, c.q. http://www.ncbi.nlm.nih.gov/sra/SRX1016818[accn]. This data set has known INDELs (http://bioinform.github.io/huref-gs/);
2. Reads with SNPs and small indels, simulated with a tool called wgsim (https://github.com/lh3/wgsim);

## 1. Illumina whole genome sequencing data set
First, download the expected results from http://bioinform.github.io/huref-gs/ in /gcc/resources/validation_data/HuRef-GS/data/results/

Next, download sequencing data:

	cd /gcc/resources/validation_data/HuRef-GS/data/original/
	wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR200/SRR2008148/SRR2008148.sra
	wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR201/SRR2015126/SRR2015126.sra

We convert these Sequence Read Archive (*.sra) files to *.fasta files, so they can be analyzed by our pipeline. For the conversion we use sra-tools.

### Download, install, use, and remove sra-tools
We temporarily download and install sra-tools in /gcc/resources/validation_data/HuRef-GS/sra-tools/ from https://github.com/ncbi/sra-tools. The README in this repo refers to the following pre-built version:

	cd /gcc/resources/validation_data/HuRef-GS/sra-tools/
	wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-centos_linux64.tar.gz
	tar -zxvf sratoolkit.2.5.2-centos_linux64.tar.gz

Convert sra to fasta (takes a while):

	cd /gcc/resources/validation_data/HuRef-GS/data
	../sra-tools/sratoolkit.2.5.2-centos_linux64/bin/fastq-dump --split-files original/SRR2008148.sra &
	../sra-tools/sratoolkit.2.5.2-centos_linux64/bin/fastq-dump --split-files original/SRR2015126.sra &

#### Split lanes

Each .sra file corresponds to one run and resulted in two .fastq files (for each end one; _1 and _2). Each .fastq file however covers two lanes. Next step is to split the reads corresponding to their lane:

	awk 'BEGIN {FS = ":"} {
	        lane=$4 
	        print $0 > "SRR2008148_L"lane"_1.fastq"
	        for (i = 1; i <= 3; i++) {
	                getline ;
	                print $0 > "SRR2008148_L"lane"_1.fastq"
	        }
	}' < SRR2008148_1.fastq

	echo "Done with SRR2008148_1.fastq"
             

	awk 'BEGIN {FS = ":"} {
	        lane=$4 
	        print $0 > "SRR2008148_L"lane"_2.fastq"
	        for (i = 1; i <= 3; i++) {
	                getline ;
	                print $0 > "SRR2008148_L"lane"_2.fastq"
	        }
	}' < SRR2008148_2.fastq

	echo "Done with SRR2008148_2.fastq"

#### Move, gzip and md5
	
	cd /gcc/resources/validation_data/HuRef-GS/data
	
First only prepare for analysis of the SRR2008148 data:
	
	gzip < SRR2008148_L1_1.fastq > /gcc/groups/gaf/tmp03/rawdata/ngs/150806_HuRef_9999_ZZZZZZZZXX/150806_HuRef_9999_ZZZZZZZZXX_L1_ZZZZZZ_1.fq.gz &
	gzip < SRR2008148_L1_2.fastq > /gcc/groups/gaf/tmp03/rawdata/ngs/150806_HuRef_9999_ZZZZZZZZXX/150806_HuRef_9999_ZZZZZZZZXX_L1_ZZZZZZ_2.fq.gz &
	gzip < SRR2008148_L2_1.fastq > /gcc/groups/gaf/tmp03/rawdata/ngs/150806_HuRef_9999_ZZZZZZZZXX/150806_HuRef_9999_ZZZZZZZZXX_L2_ZZZZZZ_1.fq.gz &
	gzip < SRR2008148_L2_2.fastq > /gcc/groups/gaf/tmp03/rawdata/ngs/150806_HuRef_9999_ZZZZZZZZXX/150806_HuRef_9999_ZZZZZZZZXX_L2_ZZZZZZ_2.fq.gz &


#### Cleaning up
As tools shouldn't pollute our resources, please delete the tool:

	cd /gcc/resources/validation_data/HuRef-GS/sra-tools
	rm -rf sratoolkit.2.5.2-centos_linux64

### Create worksheet


## 2. wgsim
Commands listed below. For details see http://biotechies.blogspot.nl/2011/12/pipelines-to-simulate-and-detect-indels.html

	cd /gcc/resources/validation_data/wgsim/ # please create if non-existent

Create small test set

	module load samtools/0.1.19
	samtools faidx /gcc/resources/b37/indices/human_g1k_v37.fasta 20 > human_g1k_v37_chr20.fasta

Index the set:

	module load bwa/0.7.12
	bwa index human_g1k_v37_chr20.fasta 

Load wgsim module

	module load wgsim/0.3.1-r13

Simulate reads:

	wgsim -N 20000000 -X 0.95 human_g1k_v37_chr20.fasta out.read1.fq out.read2.fq > wgsim.out

For illumina we expect a 1% sequencing error rate and 100 read length:
	wgsim -e 0.01 -N 20000000 -X 0.95 -1 100 -2 100 human_g1k_v37_chr20.fasta out.read1.fq out.read2.fq > wgsim.out

Improve default quality of simulated reads:

	awk ' /222/ { gsub("2", "I"); print $0; next } { print } ' out.read1.fq > out.read1b.fq
	awk ' /222/ { gsub("2", "I"); print $0; next } { print } ' out.read2.fq > out.read2b.fq 

Now align reads (According to to /gcc/tmp03/tools/NGS_DNA-2.2.1/parameters.csv we can use 8 cores):

	bwa aln -t 8 human_g1k_v37_chr20.fasta out.read1b.fq -f out1.sai
	bwa aln -t 8 human_g1k_v37_chr20.fasta out.read2b.fq -f out2.sai

Create BAM:

	bwa sampe -f out.sam -r '@RG\tID:foo\tSM:bar' human_g1k_v37_chr20.fasta out1.sai out2.sai out.read1b.fq out.read2b.fq
	samtools view -hbS out.sam -o out.bam
	samtools sort out.bam out_sorted

Call indels:

	module load delly/v0.6.7
	delly -t DEL -x human.hg19.excl.tsv -o wgsim.delly.vcf -g human_g1k_v37_chr20.fasta out_sorted.bam

Delly's output:

	0%   10   20   30   40   50   60   70   80   90   100%
	|----|----|----|----|----|----|----|----|----|----|
	***************************************************
	No structural variants found!

Manual inspection of simulated indels reveals that largest DEL was 169 bp long. Which apparently is below delly's detection limit.