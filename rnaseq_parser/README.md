

parse_rnaseq_output
--------

Works with Python 3

Python command line tool to parse all data from our PublicRNAseq pipeline into a molgenis database. Takes as input the folders with the .sh job files that were submitted to the cluster and the .out and .err log files. Need a molgenis server running locally or remotely, see https://github.com/molgenis/molgenis for more info.

For help do

```
python run_parser.py -h
```

To run the program, change the paths to the relevnat files in the CONFIG file in RNAseqParser. These options can also be changed through the command line.
