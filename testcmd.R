#Read script arguments
cargs <- commandArgs(TRUE)
args=NULL
if(length(cargs)>0){
    flags = grep("^--.*",cargs)
    values = (1:length(cargs))[-flags]
    args[values-1] = cargs[values]
    if(length(args)<tail(flags,1)){
        args[tail(flags,1)] = NA
    }
    names(args)[flags]=cargs[flags]
}

if(is.null(args) || is.na(args['--in'])) {
    print( "Usage: Rscript SNP_quality_histograms.R --nbins <number of bins (default 100)> --in <a vcf file> [--out <SNP_quality_histograms.pdf>]")
    q()
}
