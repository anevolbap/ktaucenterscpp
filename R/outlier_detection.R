outliers_tau_cutoff <- function(cutoff, mscale_bp) {
  function(ktau) {
    thr <- qchisq(cutoff, df = ktau$p)
    robust_scale <- Mscale(u = ktau$di,
                           b = mscale_bp,
                           cc = normal_consistency_constants(ktau$p))
    
    ktau$outliers <- which(ktau$di ^ 2 > thr * robust_scale ^ 2)
    
    return(ktau)
  }
  
}