###################################################################################################
# This script bins reads per chromosome in bins of  size 50.000. Input is a .bed file and output
# tsv file containing the matrix.
###################################################################################################

#Reads in commandline arguments (supplied by compute)
args<-commandArgs(TRUE)

#sets the column names for the sample for easy acces later
colnames <-c("chrom", "startpos", "stoppos", "readname", "score", "strand")
#Reads the sample, args1 = filename
sample <-(read.delim(args[1], header = FALSE, sep = "\t", quote = "\"",
                     dec = ".", fill = TRUE, col.names = colnames))
#default binsize = 50.000
bin.size = 50000
#Gets the chromosomes from all the reads from the sample as a vector. This is used to assign the read to the 
#correct chromosome
chr = as.vector(sample$chrom)
#Gets the startpositon from all the reads from the sample as a vector. This is used to assign the read to the 
#correct bin
pos = sample$startpos
#Gets a set of all the chromosomes (24 total, 1 -22 + X and Y )
chr.focus = unique(chr)[1:24]
#Removes all reads that are not on chromosome 1-22 and X or Y
index.remove = which(! chr %in% chr.focus)
chr = chr[-index.remove]
pos = pos[-index.remove]

#Gets the maximium start position of the reads, this is used to determine the number of bins
max.pos = max(pos)
#Makes an empty matrix to store the bin counts in 
bin = matrix(0, nrow = length(unique(chr)), ncol = max.pos %/% bin.size + 1, dimnames = list(unique(chr), NULL))
#Bins the reads per chromosome. Method used is integer division  by binsize + 1, for instance a read with startposition 125.000
#divided by 50.000 yields 2, + 1 is 3. The read belongs in bin 3
for (i in 1:length(pos))
{
  bin[chr[i], pos[i] %/% bin.size + 1] = bin[chr[i], pos[i] %/% bin.size + 1] + 1
}

#Writes the output, a file with the bin counts and an heatmap like image
write.table(bin, args[3], quote = FALSE, sep ="\t", row.names = TRUE)
pdf(args[2])
image(log(bin), axes = F)
axis(1,at=seq(0,1,l=length(chr.focus)),labels=chr.focus)
dev.off()
quit(save = "no", status = 0)

