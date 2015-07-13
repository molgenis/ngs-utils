#!/usr/bin/env Rscript

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

if(is.null(args) || is.na(args['--in']) || is.na(args['--step']) || is.na(args['--name'])){
	print("Usage: Rscript extract_info_GATK_variantEval.R  --in <variantEval_file> --step stepname --name dataset [--header] [--comp compRod]")
	print("Important: You need to have the GATK gsalib compiled and in your R libpath before using this script. See the GATK wiki on VariantEval to see how.")
	q()
}

library(gsalib)
d=gsa.read.gatkreport(args['--in'])

	if(!is.na(args['--header'])){
		print("Dataset,Step,#SNPs,%dbSNP,Ti/Tv Known,Ti/Tv Novel,%All comp_het_called_het,%Known comp_het_called_het,%Non-ref sensitivity,%Non-ref discrepancy,%Concordance (iChip)")
	}

	compOverlap_dbsnp = which(d$CompOverlap$CompRod=="dbSNP")
	compOverlap_all = which(d$CompOverlap$Novelty=="all")
	titv_dbsnp = which(d$TiTvVariantEvaluator$CompRod=="dbSNP")
	titv_known = which(d$TiTvVariantEvaluator$Novelty=="known")
	titv_novel = which(d$TiTvVariantEvaluator$Novelty=="novel")
	
	if(!is.na(args['--comp'])){
		geno_comp = which(d$GenotypeConcordance.simplifiedStats$CompRod==args["--comp"])
		geno_novelty = which(d$GenotypeConcordance.simplifiedStats$Novelty=="all")
		geno_row = which(d$GenotypeConcordance.simplifiedStats$row=="allSamples")
		geno_known = which(d$GenotypeConcordance.simplifiedStats$Novelty=="known")
	}
	


print(sprintf("%s,%s,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f",args['--name'],args['--step'],
d$CompOverlap$nEvalVariants[intersect(compOverlap_dbsnp,compOverlap_all)],
d$CompOverlap$compRate[intersect(compOverlap_dbsnp, compOverlap_all)],
d$TiTvVariantEvaluator$tiTvRatio[intersect(titv_dbsnp,titv_known)],
d$TiTvVariantEvaluator$tiTvRatio[intersect(titv_dbsnp,titv_novel)],
d$GenotypeConcordance.simplifiedStats$percent_comp_het_called_het[intersect(geno_comp,intersect(geno_novelty,geno_row))],
d$GenotypeConcordance.simplifiedStats$percent_comp_het_called_het[intersect(geno_comp,intersect(geno_known,geno_row))],
d$GenotypeConcordance.simplifiedStats$percent_non_reference_sensitivity[intersect(geno_comp,intersect(geno_novelty,geno_row))],
d$GenotypeConcordance.simplifiedStats$percent_non_reference_discrepancy_rate[intersect(geno_comp,intersect(geno_novelty,geno_row))],
d$GenotypeConcordance.simplifiedStats$percent_overall_genotype_concordance[intersect(geno_comp,intersect(geno_novelty,geno_row))]),
quote=FALSE)