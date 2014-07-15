args<-commandArgs(TRUE)
colnames <-c("chrom", "startpos", "stoppos", "readname", "score", "strand")


sample <-(read.delim(args[1], header = FALSE, sep = "\t", quote = "\"",
                     dec = ".", fill = TRUE, col.names = colnames))

bin.size = 50000

chr = as.vector(sample$chrom)
pos = sample$startpos

chr.focus = unique(chr)[1:24]

index.remove = which(! chr %in% chr.focus)
chr = chr[-index.remove]
pos = pos[-index.remove]

max.pos = max(pos)

bin = matrix(0, nrow = length(unique(chr)), ncol = max.pos %/% bin.size + 1, dimnames = list(unique(chr), NULL))

for (i in 1:length(pos))
{
  bin[chr[i], pos[i] %/% bin.size + 1] = bin[chr[i], pos[i] %/% bin.size + 1] + 1
}




#Writes the output
write.table(bin, args[3], quote = FALSE, sep ="\t", row.names = TRUE)



pdf(args[2])
image(log(bin), axes = F)
axis(1,at=seq(0,1,l=length(chr.focus)),labels=chr.focus)
dev.off()

quit(save = "no", status = 0)

