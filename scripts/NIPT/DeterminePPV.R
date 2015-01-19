###################################################################################################
# This script calculates the positive predictive value 


#Libs that make PPV .pdf tables
library(gridExtra)
library(lattice)
#Gets the .csv diagnostic tables
GetFileNames <- function(workdir)
{
  setwd(workdir)
  filenames <- list.files( pattern = "*DiagnostiekOutput.csv")
  return (filenames)
}
#Loads the files
LoadFiles <- function(filename)
{
  zscores <- read.delim(filename, header = TRUE, sep = ",", quote = "\"",
                        dec = ".", fill = TRUE)
  return (zscores)
}
#Function used for the integral
GetZExp <- function(x, z.obs, upper.lim, lower.lim) 
{
  return(exp(-(( x -z.obs) ^2 / 2) ) / (upper.lim - lower.lim))
}
#Function that gets the fetal upper and lower limits
GetFetalUpperLower <- function(fetal, vc)
{
  return ((fetal * 0.5) / vc)
}
CalculatePPV <- function(fetal.high, fetal.low, vc, z.obs, apriori)
{
  #Gets upper fetal limit
  upper <- GetFetalUpperLower(fetal.high, vc)
  #Gets fetal lower limit
  lower <- GetFetalUpperLower(fetal.low, vc)
  #Calculates the integral
  fetal.perc <-integrate(GetZExp, lower, upper, lower.lim = lower, upper.lim = upper , z.obs = z.obs)
  #Multiplies the integeral with the apriori risk
  fetal.apriori <- fetal.perc$value * apriori
  #Calculates PPV fractions
  ppv.frac <- fetal.apriori + (1 - apriori) * exp(-(z.obs)^2/2)
  #Converts fraction to percentage
  ppv.perc <- ((fetal.apriori / ppv.frac) * 100)
  
  return(ppv.perc)
}


GetPPV <- function(zscore, fetal.high.wide, fetal.low.wide, apriori, risk, fetal.high.narrow, fetal.low.narrow)
{
  #Logical vector to determine which PPv scores are used for median 
  normal.distributed <- vector(mode="logical", length= 4)
  #For every row in table
  for (j in 1:nrow(zscore))
  {
    #If normal distributed is yes
    if (zscore[j,2] == "Yes")
    {
      normal.distributed[j] <- TRUE
    }
    #Gets variaten coefficient
    vc <- zscore[j,1]
    #Gets observed Z score
    z.obs <- zscore[j,3]
    ppv.perc.wide <- CalculatePPV(fetal.high.wide, fetal.low.wide, vc, z.obs, apriori)
    ppv.perc.narrow <- CalculatePPV(fetal.high.narrow, fetal.low.narrow, vc, z.obs, apriori)
    ppv.perc <- (ppv.perc.wide * 0.4 + ppv.perc.narrow * 0.6)  
    #Adds PPv score to table
    zscore[j,4] <- round(ppv.perc,3) 
    names(zscore)[4] <- "PPV (%)"
  }
  #gets all ppv scores in a vector
  median.values <- zscore[,4]
  #Determines median, only normal distributed values are included
  zscore[1,5] <- median(median.values[normal.distributed])
  #Adds a priori risk to table
  zscore[1,6] <- risk
  #Decent names for new columns
  names(zscore)[5:6] <- c("Median PPV (%)"," Apriori Risk")
  #Removes NA values from table 
  zscore[is.na(zscore)] <- c("")
  return(zscore)
}
#Gets the chromosome to use as a header for the tables
GetChromosome <- function(filename)
{
  pattern <- ".+Chromosome"
  pattern2 <- gsub(pattern, "", filename)
  chromo.focus <-substr(pattern2, 0, 2)
  return(chromo.focus)
}
#This functions converts the tables to grobs 
PreparePPVTableForPrint <- function(ppvtable, chromo.focus)
{
  d <- head(ppvtable)
  table <- tableGrob(d)
  h <- grobHeight(table)
  w <- grobWidth(table)
  title <- textGrob(paste("Chromosome " , chromo.focus, sep = " "), y=unit(0.5,"npc") + 0.5*h, 
                    x=unit(0.40,"npc") + 0.5*h,
                    vjust=0, gp=gpar(fontsize=24))
  
  gt <- gTree(children=gList(table, title))
  return(gt)
}

#################################################################################################################

#Script
#Stores the command line arguments in a vector
#args[1] = temp directory where files produced during run of pipeline are stored
#args[2] = a priori risk chromosome 13 in integer
#args[3] = a priori risk chromosome 18 in integer
#args[4] = a priori risk chromosome 21 in integer
#args[5] = sample ID
args<-commandArgs(TRUE)
#Workdir variable
workdir <- args[1]
#A priori risks as integers
risk13 <- as.numeric(args[2])
risk18 <- as.numeric(args[3])
risk21 <- as.numeric(args[4])

risk.vector <- c(risk13, risk18, risk21)
#range of fetal percentages
fetal.low.wide <- 1
fetal.high.wide <- 23

fetal.low.narrow <- 6
fetal.high.narrow <- 18

#vector for apriori risks in fractions
apriori <- c( 1 / risk13, 1 / risk18, 1 / risk21) 
#Gets file names
filenames <- GetFileNames(workdir)
#Loads filenames
files <- lapply(filenames, LoadFiles)
#List to store the tables
ppv.tables <- list()
#For each table, the PPV and median PPV are calculated
for (k in 1:length(files))
{
  ppv.tables[[k]] <- GetPPV(zscore = files[[k]], fetal.high.wide = fetal.high.wide, fetal.low.wide = fetal.low.wide, apriori = apriori[k], risk = risk.vector[k],
                            fetal.high.narrow = fetal.high.narrow, fetal.low.narrow = fetal.low.narrow)
}
#List to store the tables, but now converted to grobs
ppv.pdf <- list()
#For each table, convert to grob and write .csv to disk
for (i in 1:length(ppv.tables))
{
  chromo.focus <- GetChromosome(filenames[i])
  ppv.pdf[[i]] <- PreparePPVTableForPrint(ppv.tables[[i]], chromo.focus)
  #Writes the new tables as .csv
  write.table(ppv.tables[[i]], filenames[i], quote = FALSE, sep =",", row.names = FALSE,
              col.names = TRUE)
}
#prints the pdf
pdf(paste(args[5], "_Diagnostiek_Resultaat.pdf", sep=""), height=8, width=14)
do.call(grid.arrange, c(ppv.pdf, main=args[5]))
dev.off()
