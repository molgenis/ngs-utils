#'@export
calculatePPV <- function(nipt_result, apriori)
{
  UseMethod("calculatePPV", nipt_result)
}
#'@export
calculatePPV.NCVResult <- function(nipt_result, apriori){
  print("NCVscore")
}
#'@export
calculatePPV.ZscoreResult <- function(nipt_result, apriori){
  print("Zscore")
}
#'@export
calculatePPV.RegressionResult <- function(nipt_result, apriori){
  apply(X = nipt_result$prediction_statistics, MARGIN = 2, FUN = function(x)
    GetPPV(z.obs = as.numeric(x[1]), apriori = 1/apriori, cv = as.numeric(x[2])))
}

