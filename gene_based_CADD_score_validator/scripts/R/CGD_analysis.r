library("stringr")
library("VariantAnnotation") 
library("getopt")
library("ggplot2")

#############
# multiplot #
#############
# Allows for multiplotting based upon the number of columns specified
# Expects a number of plots to draw.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) 
  {
    print(plots[[1]])
    
  } 
  else
  {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


analyse.Gene <- function(gene, vcf, plots, plotout, patient.cadd, max.thresholds){
  gene.datatable <- as.data.frame(info(vcf)[1:nrow(info(vcf)), 1:11][info(vcf)$GENE_SYMBOL == gene,])
  colours <- NULL
  results.df <- data.frame(
    "gene symbol" = gene
  )
  print(gene)
  #get non pathogenic and pathogenic indici
  index.benign = which( .02 < gene.datatable["GoNL_AF"] )
  index.benign = c(index.benign, which(.02 < gene.datatable["Thousand_Genomes_AF"]))
  index.benign = c(index.benign, which(.02 < gene.datatable["EXAC_AF"]))

  index.pathogenic = which( "5" == gene.datatable["CLINVAR_CLNSIG"] )

  # set all colour labels to unknown 
 
  print(which( "Pathogenic" %in% gene.datatable["CLINVAR_CLNSIG"] ))
  colours = rep("Unknown", nrow(gene.datatable))
  
  validGene<-F
  
  for (i in 1:nrow(gene.datatable["CLINVAR_CLNSIG"]))
  {
    elt = gene.datatable$CLINVAR_CLNSIG[[i]]
    if (0 != length(elt))
    {
      if (is.element("5", strsplit(elt, "\\|")[[1]])){
        index.pathogenic = c(index.pathogenic, i)  
      }
      else if(is.element("4", strsplit(elt, "\\|")[[1]])){
        index.pathogenic = c(index.pathogenic, i)
      }
      else if(is.element("3", strsplit(elt, "\\|")[[1]])){
        index.benign = c(index.benign, i )
      }
    }
  } 
  
  
  if(length(index.pathogenic) > 5  & length(index.benign) > 5  & (length(index.benign) + length(index.pathogenic)) > 20){
    
    
    
    # if there are any indexes found indicating pathogenic variants then set the pathogenic colour labels
    if(length(index.pathogenic) > 0)
    {
      colours[index.pathogenic] <- "Pathogenic"
    }
    # if there are any indexes found indicating Benign variants then set the Benign colour labels
    if(length(index.benign) > 0)
    {
      colours[index.benign] <- "Benign"
    }

    
    # set colour labels to the gene table
    gene.datatable$colour <- colours
    # sort data by cadd score
    gene.datatable.sorted.cadd.scaled <- gene.datatable[order(unlist(gene.datatable["CADD"])),]
    
    for(max.errors in 0:max.thresholds){
      max.benign.errors <- (length(index.benign) / 100 ) * max.errors
      max.pathogenic.errors <- (length(index.pathogenic) / 100  ) * max.errors

      cadd.scaled.benign.index <- treshold.benign.calculate(unlist(gene.datatable.sorted.cadd.scaled["colour"]), err.max=max.benign.errors)
      cadd.scaled.pathogenic.index <- treshold.pathogenic.calculate(unlist(gene.datatable.sorted.cadd.scaled["colour"]), err.max=max.pathogenic.errors)

      benign.cutoff <- unlist(gene.datatable.sorted.cadd.scaled$CADD_SCALED)[cadd.scaled.benign.index]
      pathogenic.cutoff <- unlist(gene.datatable.sorted.cadd.scaled$CADD_SCALED)[cadd.scaled.pathogenic.index]
      
      if (plots){
        
        cadd.scaled.plot <- NULL
        cadd.scaled.plot <- qplot(unlist(gene.datatable.sorted.cadd.scaled["CADD_SCALED"]), xlab="CADD scaled score", ylab="Counts", main=paste("Scaled CADD score distribution permitting", max.errors, "% errors for gene:", gene), fill=factor(unlist(gene.datatable.sorted.cadd.scaled$colour))) + scale_fill_manual(values=c("#009E73", "#D55E00", "#0072B2"))
        cadd.scaled.plot <- cadd.scaled.plot + geom_vline(xintercept = patient.cadd, colour="red", linetype = "longdash")
        
        cadd.scaled.plot <- cadd.scaled.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.cadd.scaled$CADD_SCALED)[cadd.scaled.benign.index], colour="purple")
        cadd.scaled.plot <- cadd.scaled.plot + geom_vline(xintercept = unlist(gene.datatable.sorted.cadd.scaled$CADD_SCALED)[cadd.scaled.pathogenic.index], colour="orange")
        cadd.scaled.plot <- cadd.scaled.plot + guides(fill=guide_legend(title="Classification"))
        # pdf(paste("/Users/molgenis/Desktop/afstuderen/presentatie/figures", gene, max.errors, ".pdf", sep=""),width=15, height=7)
        
        plot(cadd.scaled.plot)
        # dev.off()
      }
      
     
      if(benign.cutoff < pathogenic.cutoff){
        
        validGene = T
        results.df["Number of pathogenic variants"] = length(index.pathogenic)
        results.df["Number of benign variants"] = length(index.benign)
        
        results.df[paste("pathogenic threshold,", max.errors, "% errors")] <- pathogenic.cutoff
        results.df[paste("benign threshold,", max.errors, "% errors")] <- benign.cutoff
        
        npathogenic.pathogenic.region <- (length(which("Pathogenic" == unlist(gene.datatable.sorted.cadd.scaled$colour)[cadd.scaled.pathogenic.index: length(unlist(gene.datatable.sorted.cadd.scaled$colour))])) / (length(index.benign) + length(index.pathogenic))) * 100
        npathogenic.benign.region <- (length(which("Pathogenic" == unlist(gene.datatable.sorted.cadd.scaled$colour)[0:cadd.scaled.benign.index])) / (length(index.benign) + length(index.pathogenic))) * 100
        nbenign.pathogenic.region <- (length(which("Benign" == unlist(gene.datatable.sorted.cadd.scaled$colour)[cadd.scaled.pathogenic.index: length(unlist(gene.datatable.sorted.cadd.scaled$colour))])) / (length(index.benign) + length(index.pathogenic))) * 100
        nbenign.benign.region <-  (length(which("Benign" == unlist(gene.datatable.sorted.cadd.scaled$colour)[0:cadd.scaled.benign.index]))/ (length(index.benign) + length(index.pathogenic))) * 100
        
        results.df[paste("% pathogenic variants in pathogenic region", max.errors, "% errors")] <- npathogenic.pathogenic.region
        results.df[paste("% benign variants in pathogenic region", max.errors, "% errors")] <- nbenign.pathogenic.region
        results.df[paste("% pathogenic variants in benign region", max.errors, "% errors")] <- npathogenic.benign.region
        results.df[paste("% benign variants in benign region", max.errors, "% errors")] <- nbenign.benign.region
     }
     else{
         results.df["Number of pathogenic variants"] = length(index.pathogenic)
         results.df["Number of benign variants"] = length(index.benign)
       
         results.df[paste("pathogenic threshold,", max.errors, "% errors")] <- "N/A"
         results.df[paste("benign threshold,", max.errors, "% errors")] <- "N/A"
         
         results.df[paste("% pathogenic variants in pathogenic region", max.errors, "% errors")] <- "N/A"
         results.df[paste("% benign variants in pathogenic region", max.errors, "% errors")] <- "N/A"
         results.df[paste("% pathogenic variants in benign region", max.errors, "% errors")] <- "N/A"
         results.df[paste("% benign variants in benign region", max.errors, "% errors")] <- "N/A"
      }
    
    }
    if(validGene){
      return(results.df)
    }else{
      return(NULL)
    }
    
  }
}


treshold.benign.calculate <- function(population, err.max){
  err.counts = 0
  index = length(population)

  #check when the first uncertainty occurs
  for(i in 1:length(population)){
    if( "Pathogenic" %in% population[i]){
      err.counts <- err.counts + 1
      
      if(err.counts > err.max){
        index <- i - 1
        break
      }
      
    }
  }
 
  if ((index < 1) | (length(index) == 0)){   
    index <- 1
  }
  return(index)
}


treshold.pathogenic.calculate <- function(population, err.max){
  err.counts = 0
  index = 1

  #check when the first uncertainty occurs
  for(i in length(population):1){
    if("Benign" %in% population[i]){
      err.counts <- err.counts + 1
      
      if(err.counts > err.max){
        index <- i + 1

        break
      }
      
    }
  }

  if((index > length(population)) | (length(index) == 0)){
   
    index <- length(population)
  }

  return(index)
}


###################
# CMDline input check #
###################
spec <- matrix(c(
  'input'     , 'i', 1, "character", "input VCF file annotated using the molgenis annotation software, the following fields are required in the info collumn:        
  CADD_SCALED               Float 
  GoNL_GTC                  String
  GoNL_AF                   String
  Thousand_Genomes_AF       String
  EXAC_AF                   String
  CLINVAR_CLNSIG            String
  CLINVAR_CLNALLE           String
  GENE_SYMBOL               String",
  'out'     , 'o', 1, "character", "ouput CSV path + filename of annotated vcf-file analysis",
  'plot'    , 'p', 1, "character", "T for generating plot output and F for no plot output, defaults to false",
  'pout'    , 'q', 1, "character", "Ouput directory to write png graphs of the analysis, every gene will generate 1 plot. gene name will be used as file name",
  'help'   , 'h', 0, "logical",   "This help",
  'maxthresholds'     , 't', 1, "integer", "Defines maximum amount of unexpected values permitted to define treshold values in percentages for each class increments by 1%,
  5 will mean 0%, 1%, 2%, 3%, 4% and 5% permitted unexpected values",
  'cadd'     , 'c', 1, "double", "Defines CADD scaled score value indicating a patient"
),ncol=5,byrow=T)

# parse commandline option
opt = getopt(spec);

# check all arguments are present or help option is present, if needed show usage
if(!is.null(opt$help) || is.null(opt$input) || is.null(opt$out) || is.null(opt$cadd) || is.null(opt$maxthresholds)){
  cat(paste(getopt(spec, usage=T),"\n"));
  q()
}

# set default plot options if needed
if (is.null(opt$plot)){
  opt$plot <- F
}

# check for plot output dir
if (opt$plot){
  if(is.null(opt$pout)){
    cat(paste(getopt(spec, usage=T),"\n"));
    q()
  }
}


##########################
# annotated VCF analysis #
##########################
# start analysis
vcf.file<-readVcf(opt$input, "hg19")
genes <- unique(info(vcf.file)$GENE_SYMBOL)
cat(paste(c("## Starting analysis of: ", opt$input)))
results.rows.df<-lapply(genes, analyse.Gene, vcf=vcf.file, plots=opt$plot, plotout=opt$pout, max.thresholds=opt$maxthresholds)
df.results<-do.call("rbind", results.rows.df)
cat("## Output will be written to: ")
cat(paste(c(opt$out,"\n"), sep=""))
write.table(df.results, file = opt$out, sep = "\t", col.names = NA, qmethod = "double")
cat("#### Done ######\n")
cat(paste(c("Kept"), nrow(df.results), "out of", length(genes), "genes\n" ))

#test
vcf<-readVcf("/Users/molgenis/Desktop/Analysis/data/merged.vcf", "hg19")
#genes <- unique(unlist(info(vcf)$GENE_SYMBOL))
#length(genes)

#genes <- unique(unlist(info(vcf)$GENE_SYMBOL))
#results<-lapply(genes[2001:length(genes)], analyse.Gene, vcf=vcf, plots=F, patient.cadd=15, max.thresholds=5)
#results.df<-do.call("rbind", results)
#analyse.Gene(gene = "SMPD1", vcf = vcf, plots = T, max.thresholds = 5, patient.cadd = 10)


a <- as.data.frame(info(vcf)[1:nrow(info(vcf)), 1:11][info(vcf)$GENE_SYMBOL == "SMPD1",])

which(a$CLINVAR_CLNSIG == "2")
