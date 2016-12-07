############
# QTL_pval #
############

# function to get the QTL decomposed genetic effect for the cross-specific
# parental and ancestral linear models.


QTL_pval <- function(mppData, model, Q.eff, x, par.clu) {
  
  index <- which(substr(rownames(summary(model)$coefficients),1,3)=="QTL")
  coeffs <- summary(model)$coefficients[index, , drop = FALSE]
  
  if (Q.eff == "cr") {
    
    coeffs2 <- matrix(1, mppData$n.cr, 4)
    Q.ind <- substr(rownames(coeffs), 6, nchar(rownames(coeffs)))
    coeffs2[(unique(mppData$cross.ind) %in% Q.ind), ] <- coeffs
    
    
  } else if (Q.eff == "par") {
    
    coeffs2 <- matrix(1, mppData$n.par, 4)
    Q.ind <- substr(rownames(coeffs), 4, nchar(rownames(coeffs)))
    coeffs2[(mppData$parents %in% Q.ind), ] <- coeffs
    
    pval <- coeffs2[, 4] * sign(coeffs2[, 1])
    
    
  } else if (Q.eff=="anc") {
    
    n.anc <- length(unique(par.clu[x, ]))
    coeffs2 <- matrix(1, n.anc, 4)
    Q.ind <- substr(rownames(coeffs), 12, nchar(rownames(coeffs)))
    coeffs2[(unique(par.clu[x, ]) %in% as.numeric(Q.ind)), ] <- coeffs
    
  }
  
  pval <- coeffs2[, 4] * sign(coeffs2[, 1])
  
  if(Q.eff == "anc") {
    
    # distribute the ancestor p-value
    
    A.allele <- as.factor(par.clu[x, ])
    A <- model.matrix(~ A.allele - 1)
    pval <- as.vector(A %*% pval)
    
  }
  
  return(pval)
  
}