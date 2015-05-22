# Steps to auto-validate our Variant Calling pipeline by spiking-in simulated reads
Rough idea is to simulate reads from phiX with known SNPs, spike-in those reads at the start of our pipeline, and if successful our pipeline should find and QC-report those SNPs.

## 1: sim software
We have chosen to use 'wgsim' (`https://github.com/lh3/wgsim`) to simulate reads.

Download the software
```bash
	wget https://github.com/lh3/wgsim/archive/master.zip
```
unzip it and compile
```bash
	gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
```