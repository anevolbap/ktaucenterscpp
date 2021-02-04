MSCALE_BREAKDOWN_POINT = 0.5

flag_outliers <- function(ktau, cutoff) {
  
  thr <- qchisq(cutoff, df = ktau$p)
    
  robust_scale <- Mscale(u = ktau$di,
                         b = MSCALE_BREAKDOWN_POINT,
                         cc = normal_consistency_constants(ktau$p))
  
  outliers <- which(ktau$di^2 > thr * robust_scale^2)
  
  return(outliers)
}