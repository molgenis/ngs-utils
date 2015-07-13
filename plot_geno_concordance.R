
args <- commandArgs(TRUE)
if(is.na(args[1]) || is.na(args[2])){
	print( "Usage: Rscript plot_geno_concordance.R <concordance_file> <out_plot.jpg> [plot_title]")
	q()
}
title = "Data concordance"
if(!is.na(args[3])){
	title = args[3]
}

con = read.csv(args[1], sep="\t", header=TRUE)
con$concordance = 100-con$discordant...
con.ordered = con[order(con$concordance),]

jpeg(args[2], width=700, height=500)
barplot(con.ordered$concordance, names=con.ordered$file, las=2, main=sprintf("%s\nmean=%.3f, med=%.3f",title, mean(con$concordance), median(con$concordance)), ylab=sprintf("Concordance (%%) over ~%dK loci from Immunochip",round(median(con$loci.in.both.files)/1000)), ylim=c(96,100), xpd=FALSE)
abline(h = mean(con$concordance), col="orange")
dev.off()

