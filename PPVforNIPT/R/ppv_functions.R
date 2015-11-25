GetZExp <- function(x, z.obs, upper.lim, lower.lim) {
  return(exp(-(( x -z.obs) ^2 / 2) ) / (upper.lim - lower.lim))
}
GetFetalUpperLower <- function(fetal, vc){
  return ((fetal * 0.5) / vc)
}
CalculatePPV <- function(fetal.high, fetal.low, vc, z.obs, apriori){
  upper <- GetFetalUpperLower(fetal.high, (vc *100))
  lower <- GetFetalUpperLower(fetal.low, (vc*100))
  fetal.perc <-integrate(GetZExp, lower, upper, lower.lim = lower, upper.lim = upper , z.obs = z.obs)
  fetal.apriori <- fetal.perc$value * apriori
  ppv.frac <- fetal.apriori + (1 - apriori) * exp(-(z.obs)^2/2)
  ppv.perc <- ((fetal.apriori / ppv.frac) * 100)
  return(round(ppv.perc, digits = 3))
}
GetPPV <- function(z.obs, apriori, cv){
  ppv.perc.wide <- CalculatePPV(fetal.high.wide, fetal.low.wide, cv, z.obs, apriori)
  ppv.perc.narrow <- CalculatePPV(fetal.high.narrow, fetal.low.narrow, cv, z.obs, apriori)
  ppv.perc <- (ppv.perc.wide * 0.4 + ppv.perc.narrow * 0.6)

  return(ppv.perc)
}
