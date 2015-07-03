
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
	print("Usage: Rscript plot_shared_loci.R  --in <concordance_file> --out <out_plot.jpg> [options]")
	print("Options:")
	print("--data1 NAME: Name of the first dataset.")
	print("--data2 NAME: Name of the second dataset")
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

a = read.csv(args['--in'], sep="\t", header=TRUE)
b = a[order(a$loci.in.both.files),]
c = data.frame(b$loci.in.both.files,b$loci.in.file.2.only,b$loci.in.file.1.only)
d =t(as.matrix(c))
ord = colSums(d)
jpeg(args['--out'], height=500, width=700)
par(xpd=T, mar=par()$mar+c(0,0,6,0))
barplot(d[,order(colSums(d))],names.arg=b$file,main=sprintf("%s / %s shared loci",d1,d2), ylab="Loci", col=rainbow(3),las=2)
legend("top",inset=-.185,c("Shared", sprintf("%s only",d1),sprintf("%s only",d2)), cex=0.8, fill=rainbow(3));
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()




