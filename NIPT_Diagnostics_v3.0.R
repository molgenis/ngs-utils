library(NIPTeR)
library(PPVforNIPT)
library(rmarkdown)
library(knitr)

args <- commandArgs(TRUE)
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

outputpath <- args[1]
sample_name <- args[2]
control_group_rds <- args[3]
raw_sample_bam <- args[4]
Rmarkdownfile_Diagnostic_output_table <- args[5]
Rmarkdownfile_QC <- args[6]

a_priori_13 <- args[7]
a_priori_18 <- args[8]
a_priori_21 <- args[9]

message ("Analyzing sample ", sample_name)

#read control group
message ( "Loading control samples" ) 
NIPT_control_group_separate <- readRDS( file = control_group_rds )

#Create bins in sample
NIPT_raw_sample_bin <- bin_bam_sample( bam_filepath = raw_sample_bam, separate_strands = T )

#Analysis
message ("Analyzing sample")

#Do chi squared correction on sample and control group
NIPT_chi_corrected_data <- chi_correct( nipt_sample = NIPT_raw_sample_bin, nipt_control_group = NIPT_control_group_separate )
NIPT_chi_corrected_sample <- NIPT_chi_corrected_data$sample
NIPT_chi_corrected_controls <- NIPT_chi_corrected_data$control_group

#Calculate regression-based Z-score
regression_chr_13 <- perform_regression( nipt_sample = NIPT_chi_corrected_sample, nipt_control_group = NIPT_chi_corrected_controls, chromo_focus = 13 )
regression_chr_18 <- perform_regression( nipt_sample = NIPT_chi_corrected_sample, nipt_control_group = NIPT_chi_corrected_controls, chromo_focus = 18 )
regression_chr_21 <- perform_regression( nipt_sample = NIPT_chi_corrected_sample, nipt_control_group = NIPT_chi_corrected_controls, chromo_focus = 21 )

#Retrieve coefficient of variation from regression
CV_13_1 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[2,1]))*100
CV_13_2 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[2,2]))*100
CV_13_3 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[2,3]))*100
CV_13_4 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[2,4]))*100
CV_18_1 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[2,1]))*100
CV_18_2 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[2,2]))*100
CV_18_3 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[2,3]))*100
CV_18_4 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[2,4]))*100
CV_21_1 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[2,1]))*100
CV_21_2 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[2,2]))*100
CV_21_3 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[2,3]))*100
CV_21_4 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[2,4]))*100

#Check normal distribution
Normal_13_set1 <- ifelse(as.vector(regression_chr_13$prediction_statistics[4,1])>0.05, "Yes", "No")
Normal_13_set2 <- ifelse(as.vector(regression_chr_13$prediction_statistics[4,2])>0.05, "Yes", "No")
Normal_13_set3 <- ifelse(as.vector(regression_chr_13$prediction_statistics[4,3])>0.05, "Yes", "No")
Normal_13_set4 <- ifelse(as.vector(regression_chr_13$prediction_statistics[4,4])>0.05, "Yes", "No")
Normal_18_set1 <- ifelse(as.vector(regression_chr_18$prediction_statistics[4,1])>0.05, "Yes", "No")
Normal_18_set2 <- ifelse(as.vector(regression_chr_18$prediction_statistics[4,2])>0.05, "Yes", "No")
Normal_18_set3 <- ifelse(as.vector(regression_chr_18$prediction_statistics[4,3])>0.05, "Yes", "No")
Normal_18_set4 <- ifelse(as.vector(regression_chr_18$prediction_statistics[4,4])>0.05, "Yes", "No")
Normal_21_set1 <- ifelse(as.vector(regression_chr_21$prediction_statistics[4,1])>0.05, "Yes", "No")
Normal_21_set2 <- ifelse(as.vector(regression_chr_21$prediction_statistics[4,2])>0.05, "Yes", "No")
Normal_21_set3 <- ifelse(as.vector(regression_chr_21$prediction_statistics[4,3])>0.05, "Yes", "No")
Normal_21_set4 <- ifelse(as.vector(regression_chr_21$prediction_statistics[4,4])>0.05, "Yes", "No")

#Retrieve Z-scores
RegZscore_13_1 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[1,1]))
RegZscore_13_2 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[1,2]))
RegZscore_13_3 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[1,3]))
RegZscore_13_4 <- as.numeric(as.vector(regression_chr_13$prediction_statistics[1,4]))
RegZscore_18_1 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[1,1]))
RegZscore_18_2 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[1,2]))
RegZscore_18_3 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[1,3]))
RegZscore_18_4 <- as.numeric(as.vector(regression_chr_18$prediction_statistics[1,4]))
RegZscore_21_1 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[1,1]))
RegZscore_21_2 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[1,2]))
RegZscore_21_3 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[1,3]))
RegZscore_21_4 <- as.numeric(as.vector(regression_chr_21$prediction_statistics[1,4]))

#Calculate Posterior Risk (%) 
posterior_risk_13 <- PPVforNIPT::calculatePPV(regression_chr_13, a_priori_13)
posterior_risk_18 <- PPVforNIPT::calculatePPV(regression_chr_18, a_priori_18)
posterior_risk_21 <- PPVforNIPT::calculatePPV(regression_chr_21, a_priori_21)

posterior_risk_13_1 <- posterior_risk_13[1]
posterior_risk_13_2 <- posterior_risk_13[2]
posterior_risk_13_3 <- posterior_risk_13[3]
posterior_risk_13_4 <- posterior_risk_13[4]
posterior_risk_18_1 <- posterior_risk_18[1]
posterior_risk_18_2 <- posterior_risk_18[2]
posterior_risk_18_3 <- posterior_risk_18[3]
posterior_risk_18_4 <- posterior_risk_18[4]
posterior_risk_21_1 <- posterior_risk_21[1]
posterior_risk_21_2 <- posterior_risk_21[2]
posterior_risk_21_3 <- posterior_risk_21[3]
posterior_risk_21_4 <- posterior_risk_21[4]

#Calculate median Posterior risk (%)
normal_distributed13 <- c(Normal_13_set1, Normal_13_set2, Normal_13_set3, Normal_13_set4)
normal_distributed18 <- c(Normal_18_set1, Normal_18_set2, Normal_18_set3, Normal_18_set4)
normal_distributed21 <- c(Normal_21_set1, Normal_21_set2, Normal_21_set3, Normal_21_set4)

datafr_for_median13 <- data.frame(normal_distributed13, posterior_risk_13)
datafr_for_median18 <- data.frame(normal_distributed18, posterior_risk_18)
datafr_for_median21 <- data.frame(normal_distributed21, posterior_risk_21)

datafr_for_median13_selected <- datafr_for_median13[datafr_for_median13$normal_distributed13=="Yes",]$posterior_risk_13
datafr_for_median18_selected <- datafr_for_median18[datafr_for_median18$normal_distributed18=="Yes",]$posterior_risk_18
datafr_for_median21_selected <- datafr_for_median21[datafr_for_median21$normal_distributed21=="Yes",]$posterior_risk_21

median_posterior_risk_13 <- median(datafr_for_median13_selected)
median_posterior_risk_18 <- median(datafr_for_median18_selected)
median_posterior_risk_21 <- median(datafr_for_median21_selected)

#Extra Quality Control Metrics
message( "Calculating extra QC Metrics" )

#Retrieve prediction statistics
Pred_stat_chr13 <- regression_chr_13$prediction_statistics
Pred_stat_chr18 <- regression_chr_18$prediction_statistics
Pred_stat_chr21 <- regression_chr_21$prediction_statistics

#Matchscore calculation
Match_score_list <- match_control_group( NIPT_raw_sample_bin, NIPT_control_group_separate, mode = "report" )
Average_match_score <- sum(Match_score_list) / length(Match_score_list)

#Alternative Z-score calculation (after GC correction)
gc_corrected_sample <- gc_correct(nipt_object = NIPT_raw_sample_bin, method = "bin")
gc_correctied_control_group_separate <- gc_correct(nipt_object = NIPT_control_group_separate, method = "bin" )
Zscore13 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 13 )
Zscore18 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 18 )
Zscore21 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 21 )
CV_from_Zscore13 <-unname(Zscore13$control_group_statistics[2]/Zscore13$control_group_statistics[1]*100)
CV_from_Zscore18 <-unname(Zscore18$control_group_statistics[2]/Zscore18$control_group_statistics[1]*100)
CV_from_Zscore21 <-unname(Zscore21$control_group_statistics[2]/Zscore21$control_group_statistics[1]*100)
Normal_13_fromZscore <- ifelse(as.vector(Zscore13$control_group_statistics[3])>0.05, "Yes", "No")
Normal_18_fromZscore <- ifelse(as.vector(Zscore18$control_group_statistics[3])>0.05, "Yes", "No")
Normal_21_fromZscore <- ifelse(as.vector(Zscore21$control_group_statistics[3])>0.05, "Yes", "No")

#Check if Z-scores of all chromosomes except 13, 18 and 21 are within -3 +3 boundaries
Zscore1 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 1 )
Zscore2 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 2 )
Zscore3 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 3 )
Zscore4 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 4 )
Zscore5 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 5 )
Zscore6 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 6 )
Zscore7 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 7 )
Zscore8 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 8 )
Zscore9 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 9 )
Zscore10 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 10 )
Zscore11 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 11 )
Zscore12 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 12 )
Zscore14 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 14 )
Zscore15 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 15 )
Zscore16 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 16 )
Zscore17 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 17 )
Zscore19 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 19 )
Zscore20 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 20 )
Zscore22 <- calculate_z_score(nipt_sample = gc_corrected_sample, nipt_control_group = gc_correctied_control_group_separate, chromo_focus = 22 )

All_chromosomes_normal_Zscore <- ifelse(Zscore1$sample_Zscore<3 & Zscore2$sample_Zscore<3 & Zscore3$sample_Zscore<3 & Zscore4$sample_Zscore<3 & Zscore5$sample_Zscore<3 & 
                                          Zscore6$sample_Zscore<3 & Zscore7$sample_Zscore<3 & Zscore8$sample_Zscore<3 & Zscore9$sample_Zscore<3 & Zscore10$sample_Zscore<3 & 
                                          Zscore11$sample_Zscore<3 & Zscore12$sample_Zscore<3 & Zscore14$sample_Zscore<3 & Zscore15$sample_Zscore<3 & Zscore16$sample_Zscore<3 & 
                                          Zscore17$sample_Zscore<3 & Zscore19$sample_Zscore<3 & Zscore20$sample_Zscore<3 & Zscore22$sample_Zscore<3 & 
                                          Zscore1$sample_Zscore>-3 & Zscore2$sample_Zscore>-3 & Zscore3$sample_Zscore>-3 & Zscore4$sample_Zscore>-3 & Zscore5$sample_Zscore>-3 & 
                                          Zscore6$sample_Zscore>-3 & Zscore7$sample_Zscore>-3 & Zscore8$sample_Zscore>-3 & Zscore9$sample_Zscore>-3 & Zscore10$sample_Zscore>-3 & 
                                          Zscore11$sample_Zscore>-3 & Zscore12$sample_Zscore>-3 & Zscore14$sample_Zscore>-3 & Zscore15$sample_Zscore>-3 & Zscore16$sample_Zscore>-3 & 
                                          Zscore17$sample_Zscore>-3 & Zscore19$sample_Zscore>-3 & Zscore20$sample_Zscore>-3 & Zscore22$sample_Zscore>-3,  "Yes", "No" )

#Create output tables 
message( "Creating tables" )

#Standard diagnostics output table
   rmarkdown::render(Rmarkdownfile_Diagnostic_output_table, params = list(
  samplename = sample_name, CV13_1 = CV_13_1, CV13_2 = CV_13_2, CV13_3 = CV_13_3, CV13_4 = CV_13_4, 
  Normal13_1 = Normal_13_set1, Normal13_2 = Normal_13_set2, Normal13_3 = Normal_13_set3, Normal13_4 = Normal_13_set4,
  Zscore13_1 = RegZscore_13_1, Zscore13_2 = RegZscore_13_2, Zscore13_3 = RegZscore_13_3, Zscore13_4 = RegZscore_13_4,
  posteriorrisk13_1 = posterior_risk_13_1, posteriorrisk13_2 = posterior_risk_13_2, posteriorrisk13_3 = posterior_risk_13_3, posteriorrisk13_4 = posterior_risk_13_4,
  medianposteriorrisk13 = median_posterior_risk_13, 
  apriori13 = a_priori_13,
  CV18_1 = CV_18_1, CV18_2 = CV_18_2, CV18_3 = CV_18_3, CV18_4 = CV_18_4, 
  Normal18_1 = Normal_18_set1, Normal18_2 = Normal_18_set2, Normal18_3 = Normal_18_set3, Normal18_4 = Normal_18_set4,
  Zscore18_1 = RegZscore_18_1, Zscore18_2 = RegZscore_18_2, Zscore18_3 = RegZscore_18_3, Zscore18_4 = RegZscore_18_4,
  posteriorrisk18_1 = posterior_risk_18_1, posteriorrisk18_2 = posterior_risk_18_2, posteriorrisk18_3 = posterior_risk_18_3, posteriorrisk18_4 = posterior_risk_18_4,
  medianposteriorrisk18 = median_posterior_risk_18, 
  apriori18 = a_priori_18,
  CV21_1 = CV_21_1, CV21_2 = CV_21_2, CV21_3 = CV_21_3, CV21_4 = CV_21_4, 
  Normal21_1 = Normal_21_set1, Normal21_2 = Normal_21_set2, Normal21_3 = Normal_21_set3, Normal21_4 = Normal_21_set4,
  Zscore21_1 = RegZscore_21_1, Zscore21_2 = RegZscore_21_2, Zscore21_3 = RegZscore_21_3, Zscore21_4 = RegZscore_21_4,
  posteriorrisk21_1 = posterior_risk_21_1, posteriorrisk21_2 = posterior_risk_21_2, posteriorrisk21_3 = posterior_risk_21_3, posteriorrisk21_4 = posterior_risk_21_4,
  medianposteriorrisk21 = median_posterior_risk_21, 
  apriori21 = a_priori_21
))


#QC output table
   rmarkdown::render(Rmarkdownfile_QC, params = list(
     samplename = sample_name, Z_score13 = Zscore13$sample_Zscore, Z_score18 = Zscore18$sample_Zscore, Z_score21 = Zscore21$sample_Zscore, 
     CV13 = CV_from_Zscore13, CV18 = CV_from_Zscore18, CV21 = CV_from_Zscore21,
     Normal13 = Normal_13_fromZscore, Normal18 =  Normal_18_fromZscore, Normal21 = Normal_21_fromZscore, 
     All_chromosomes_norm_Zscore = All_chromosomes_normal_Zscore, AvMatchScore = Average_match_score, MatchScorelist = Match_score_list,
     PredictionStatistics13 = Pred_stat_chr13, PredictionStatistics18 = Pred_stat_chr18, PredictionStatistics21 = Pred_stat_chr21))
