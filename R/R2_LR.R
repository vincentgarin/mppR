#########
# R2_LR #
#########

# function to compute the R squared of a list of QTL using Likelihood R squared
# using formula from Sun et al. 2010

# argument

# mppData mppData

# QTL list of QTL

# adjust should the R squared value be adjusted or not

R2_LR <- function(mppData, QTL){
  
  trait <- mppData$trait[, 1]
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  cross.ind <- mppData$cross.ind
  
  # remove non complete observations
  
  dataset <- cbind(mppData$trait[, 1], cross.mat, QTL, cross.ind)
  index <- complete.cases(dataset)
  trait <- trait[index]
  cross.mat <- cross.mat[index, , drop = FALSE] 
  cross.ind <- cross.ind[index]
  QTL <- QTL[index, , drop = FALSE]
  n <- length(trait)
  
  
  # full model
  
  model.full <- tryCatch(gls(model = trait ~ -1 + cross.mat + QTL,
                             method = "ML",
                             weights = varIdent(form = ~ 1 | cross.ind)),
                         error=function(e) NULL)
  
  # reduce model
  
  model.red <- tryCatch(gls(model = trait ~ -1 + cross.mat, method = "ML",
                            weights = varIdent(form = ~ 1 | cross.ind)),
                        error=function(e) NULL)
  
  if(is.null(model.full) || is.null(model.red)){ R2 <- R2.adj <-  NA
  
  } else {
    
    R2 <- 1-exp(-(2/n)*(model.full$logLik[1] - model.red$logLik[1]))
    
    # adjust R2
      
    z <- dim(QTL)[2]
    N <- as.numeric(substring(attr(anova(model.full), "label"), 12, 14))
    R2.adj <- R2 - ((z/N)*(1-R2))
    
    R2 <- 100 * R2; R2.adj <- 100 * R2.adj
      
    
  }
  
  return(c(R2, R2.adj))
  
}