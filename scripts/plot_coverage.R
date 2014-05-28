
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
	print("--png: outputs a PNG file rather than a PDF (Needs X11 capabilities)")
	print("NOTE: You need the R.utils package to run this script")
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

title = "Coverage"
if(!is.na(args['--title'])){
	title = args['--title']
}


#Open file and prepare data
in_file <- read.csv(args['--in'],sep="\t", row.names=1, header=file_header)
a = as.data.frame(t(as.matrix(in_file)))

#set max_depth
max_depth = dim(a)[1];
if(!is.na(args['--max-depth']) && as.integer(args['--max-depth'])>1 ){
	max_depth = as.integer(args['--max-depth'])
}

#Normalize data
a = a / a[1,1]
for (i in seq(1,max_depth+1)){ a[i,]=a[i,]-a[i+1,]}
row.names(a) <-  seq(0,length(row.names(a))-1)

#Create poisson
p = dpois(seq(0,max_depth),expected_coverage)
#p = dnorm(seq(0,max_depth),mean=expected_coverage,sd=7)
max_poiss = 0.2
if(expected_coverage > 0){
	max_poiss = ceiling(max(p)*100)/100
}else{
	p[1] = 0
}

#Plot poisson
if(is.element('--png',names(args))){
	png(args['--out'], width=800, height=800)
}else{
	pdf(args['--out'])
}
plot(p, type="l", col="blue", xlim=c(1,max_depth+1),xaxt='n',xaxs='i',yaxs='i',ylim=c(0,max_poiss), main=title, xlab="Coverage", ylab="");
axis(1,at=1:(max_depth+1),lab=seq(0,max_depth))

#Plot all curves with sufficient coverage and tag outliers
count = 1 
outliers = c()
for ( c in a ){
	if(expected_coverage > 0 && (which.max(c[2:max_depth])<expected_coverage-1)){
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
		for(i in seq(1,length(outliers))){
 			points( a[,outliers[i]], type="l",xaxs='i', col=rainb[i],lwd=2 )
			text(1.5,0.115-i*0.004,paste(outliers[i]," (max=", which.max(a[,outliers[i]][2:max_depth])-1 ,")",sep=""),col=rainb[i],cex=1,adj=c(0,0))
		}
	}
	
	#Plot Poisson again (on top)
	pcol = rainb[length(outliers)+1]
	lines(p,col=pcol,lwd=2)
	text(1.5,0.115,sprintf("Poisson (lambda=%d)",expected_coverage),col=pcol,cex=1,adj=c(0,0))
	abline(v=expected_coverage+1,col=pcol)
}

dev.off()

