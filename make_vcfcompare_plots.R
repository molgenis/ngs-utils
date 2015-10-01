library(stringr)

args=(commandArgs(TRUE))

stats.file = t(c(args[1], args[2],args[3]))

for (j in 1:nrow(stats.file))
{
    DELLYTYPE = stats.file[j,1]
    lines = readLines(stats.file[j,2])
    COMPARED = stats.file[j,3]


    xLabels = NULL
    FDR = NULL
    TPR = NULL
    for (i in 1:length(lines))
    {
        if (DELLYTYPE == lines[i]) break
    }

    while (">" != str_sub(lines[i], 2, 2)) i = i + 1
    i = i + 1

    while ("[" == str_sub(lines[i], 1, 1))
    {
        xLabels = c(xLabels, lines[i])
        vals    = str_split(lines[i], ",")[[1]]
        TPR     = c(TPR, as.numeric(str_split(vals[2], ":")[[1]][2]))
        FDR     = c(FDR, as.numeric(vals[3]))

        i = i + 1
    }

    pdf(paste0(stats.file[j,2],".",args[1],".pdf"))
        par(mai = c(4.02, 0.82, 0.82, 0.42))
        barplot(rbind(TPR, FDR), col=c('lightblue', 'mistyrose'), beside=T, legend.text=c('TPR', 'FDR'), args.legend=list(x="topright"), main=paste0("Discovery rate for: ", DELLYTYPE," ", COMPARED), xlab = "", ylab="Percentage (%)", names.arg = 
xLabels, las=2)
        title(xlab = "[length]: TPR, FDR, TP, FP, T", line=18)
    dev.off()
}
