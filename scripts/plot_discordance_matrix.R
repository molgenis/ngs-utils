

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

if(is.null(args) || is.na(args['--in']) || is.na(args['--out'])){
	print( "Usage: Rscript plot_discordance_matrix.R --in <discordance_matrix_file> --out <out_plot.jpg> [options]")
	print("Options:")
	print("--data1 NAME: Name of the first dataset.")
	print("--data2 NAME: Name of the second dataset.")
	print("--show-concordant: Show concordant loci as well.")
	print("--show-unknowns <concordance_file>: If present, the concordance file has to be over the same data to plot as 'unknown' all loci that were not captured by the concordance matrix since the alleles were not exact matches (e.g. if one of the allele was monomorphic in one set).")
	q()
}
d1 = "Dataset1"
if(!is.na(args['--data1'])){
	d1 = args['--data1']
}
d2 = "Dataset2"
if(!is.na(args['--data2'])){
	d2 = args['--data2']
}

a = read.csv(args['--in'], sep="\t",header=TRUE)
legend_texts = c(sprintf("%s 0/0 - %s 0/0",d1,d2),sprintf("%s 0/1 - %s 0/0",d1,d2),sprintf("%s 1/1 - %s 0/0",d1,d2),sprintf("%s ./. - %s 0/0",d1,d2),sprintf("%s 0/0 - %s 0/1",d1,d2),sprintf("%s 0/1 - %s 0/1",d1,d2),sprintf("%s 1/1 - %s 0/1",d1,d2),sprintf("%s ./. - %s 0/1",d1,d2),sprintf("%s 0/0 - %s 1/1",d1,d2),sprintf("%s 0/1 - %s 1/1",d1,d2),sprintf("%s 1/1 - %s 1/1",d1,d2),sprintf("%s ./. - %s 1/1",d1,d2), sprintf("%s 0/0 - %s ./.",d1,d2), sprintf("%s 0/1 - %s ./.",d1,d2),sprintf("%s 1/1 - %s ./.",d1,d2), sprintf("%s ./. - %s ./.",d1,d2))
rows=seq(2,17,1)
concordant_rows = c(2,7,12)
#If a concordance file has been passed, get the unknown factor
if(!is.na(args['--show-unknowns'])){
	x = read.csv(args['--show-unknowns'], sep="\t",header=TRUE)
	a$total =  x$discordant - rowSums(a[,-c(1,concordant_rows)])
	legend_texts = c(legend_texts,"unknown")
	rows = c(rows,18)
}

#Hide concordance unless show concordance option is selected 
if(!is.element('--show-concordant',names(args))){
	legend_texts = legend_texts[-concordant_rows+1]	
	rows = rows[-concordant_rows+1]
}
b = t(as.matrix(a[rows]))
#b = b[rows,]
colnames(b) = a$file
e = order(colSums(b))
jpeg(args['--out'], height=500, width=700)
mycol = rainbow(length(rows))
barplot(b[,e], las=2, col=mycol,main=sprintf("%s / %s genotype discordance",d1,d2), ylab="#Loci", las=2 )
abline(h=median(colSums(b)), col="orange", lwd=2)

legend("top",legend_texts, ncol=3,cex=0.8, fill=mycol);
dev.off()

