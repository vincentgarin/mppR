##########
# R2_lin #
##########

# function to compute the R squared of a list of QTL

# argument

# mppData mppData

# QTL list of QTL

# adjust should the R squared value be adjusted or not

R2_lin <- function(mppData, QTL){
  
  trait <- mppData$trait[, 1]
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  # remove non complete observations
  
  dataset <- cbind(mppData$trait, cross.mat, QTL)
  index <- complete.cases(dataset)
  trait <- trait[index]
  cross.mat <- cross.mat[index, , drop = FALSE]
  QTL <- QTL[index, , drop = FALSE]
  
  
  model.full <- tryCatch(lm(trait~-1 + cross.mat + QTL),
                         error=function(e) NULL)
  
  model.fam <- tryCatch(lm(trait~-1+cross.mat),
                        error=function(e) NULL)
  
  if(is.null(model.full) || is.null(model.fam)){R2 <- R2.adj <- NA
  
  } else {
    
    RSS.full <-  anova(model.full)$S[3]
    RSS.fam <- anova(model.fam)$S[2]
    
    R2 <- 1-(RSS.full/RSS.fam)
    
    # adjust
    
    z <- dim(QTL)[2]
    N <- anova(model.full)[3, 1]
    R2.adj <- R2 - ((z/N)*(1-R2))
    
    R2 <- 100 * R2; R2.adj <- 100 * R2.adj
    
  }
  
  return(c(R2, R2.adj))
  
}