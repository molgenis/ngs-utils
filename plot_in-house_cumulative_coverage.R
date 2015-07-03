
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
	print( "Usage: Rscript plot_coverage.R --in <coverage_file> --out <out_plot.pdf> [options]")
	print( "Options:")
	print("--max-depth DEPTH: Maximum depth of coverage to plot. If <1, then all dephts are plotted.)")
	print("--expected-coverage COVERAGE: The expected coverage. If set, the expected coverage is plotted and curves below the expected coverage are highlighted.")
	print("--no-header: The file contains no header.")
	print("--title TITLE: Sets the plot title to the given title")
	print("--png: Outputs a PNG rather than a PDF. Needs X11 capabilities.")
	q()
}

#Read arguments
expected_coverage = 0
if(!is.na(args['--expected-coverage'])){
	expected_coverage = as.integer(args['--expected-coverage'])
}

file_header = TRUE
if(is.element('--no-header',names(args))){
	file_header = FALSE
}

title = "Cumulative Coverage"
if(!is.na(args['--title'])){
	title = args['--title']
}


#Open file and prepare data
in_file <- read.csv(args['--in'],sep="\t", row.names=1, header=file_header)
a = as.data.frame(t(as.matrix(in_file)))
a = a / a[1,1]
row.names(a) <-  seq(0,length(row.names(a))-1)

#set max_depth
max_depth = dim(a)[1];
if(!is.na(args['--max-depth']) && as.integer(args['--max-depth'])>1 ){
	max_depth = as.integer(args['--max-depth'])
}

#Create and plot poisson
p = 1-ppois(seq(0,max_depth),expected_coverage)
if(is.element('--png',names(args))){
	png(args['--out'], width=800, height=800)
}else{
	pdf(args['--out'])
}
plot(p, type="l", col="blue", xlim=c(1,max_depth+1), ylim=c(0,1), xaxt='n',xaxs='i',yaxs='i',main=title, xlab="Coverage", ylab="",yaxt='n');
axis(1,at=1:(max_depth+1),lab=seq(0,max_depth))
axis(2,seq(0,1,0.1))

#Plot all curves with sufficient coverage and tag outliers
count = 1 
outliers = c()
for ( c in a ){ 
	if(expected_coverage>0 && c[expected_coverage+1]<p[expected_coverage+1]){
		outliers = c(outliers,names(a)[count])
	}
	else{
		points( c, type="l",xaxs='i' )
	}
	count=count+1
}

#Plot additionnal information if expected coverage is set
if(expected_coverage>0){
	rainb = sample(rainbow(length(outliers)+1))
	#Plot outliers
	if(length(outliers)>0){
		for (i in seq(1,length(outliers))){
			outlier =  a[,outliers[i]]
		 	points(outlier, type="l",xaxs='i', col=rainb[i],lwd=2 )
			text(expected_coverage+1.5,outlier[expected_coverage+1],round(outlier[expected_coverage+1],2),cex=1,col=rainb[i], adj=c(0,0))
			text(1.5,0.025+i*0.03,outliers[i],col=rainb[i],cex=1,adj=c(0,0))
		}
	}
	#Plot poisson
	pcol = rainb[length(outliers)+1]
	lines(p,col=pcol,lwd=2)
	text(expected_coverage+1.5,p[expected_coverage+1],round(p[expected_coverage+1],2),cex=1,col=pcol,adj=c(0,0))
	text(1.5,0.025,sprintf("Poisson (lambda=%d)",expected_coverage),col=pcol,cex=1,adj=c(0,0))
	abline(v=expected_coverage+1,col=pcol,lwd=2)
}

# plot points for minimal desired coverage
points(c(1, 10, 20), c(.92, .9, .8), pch=19, cex=1, col="darkred")
legend("topright", pch=19, legend="desired coverage", col="darkred", box.lty=0)
dev.off()

