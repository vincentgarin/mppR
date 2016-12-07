################
# QTL_pval_mix #
################

# function to get the QTL decomposed genetic effect for the cross-specific
# parental and ancestral mixed models.


QTL_pval_mix <- function(model, Q.eff, QTL.el, x, par.clu, fct) {
  
  if(fct == "SIM"){
    start.ind <- 2; end.ind <- 1
  } else if(fct == "CIM") {
    start.ind <- 3; end.ind <- 2}
  
  sign <- sign(rev(model$coefficients$fixed[1:QTL.el]))
  pval <- wald(model)[start.ind:(QTL.el + end.ind), 4]
  pval <- pval * sign
  pval[pval == 0] <- 1
  
  if (Q.eff == "par") {
    
    pval <- c(1, pval) # put a 1 for the parent that was in the intercept.
    
    
  } else if (Q.eff == "anc") {
    
    pval <- c(1, pval)  
    
    # distribute the ancestor p-value
    
    A.allele <- as.factor(par.clu[x, ])
    A <- model.matrix(~ A.allele - 1)
    pval <- as.vector(A %*% pval)
    
  }
  
  return(pval)
  
  
} 