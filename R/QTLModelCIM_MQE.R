###################
# QTLModelCIM_MQE #
###################

QTLModelCIM_MQE <- function(x, mppData, mppData_bi, cross.mat, par.mat,
                            Qeff.part, par.clu, VCOV, cof.list, cof.part){
  
  # 1. formation of the QTL incidence matrix
  ###########################################
  
  ### 2.2 QTL position
  
  if(Qeff.part[x] == "biall"){
    
    QTL <- IncMat_QTL(x = x, mppData = mppData_bi, cross.mat = cross.mat,
                      par.mat = par.mat, par.clu = par.clu, Q.eff = Qeff.part[x])
    
  } else {
    
    QTL <- IncMat_QTL(x = x, mppData = mppData, cross.mat = cross.mat,
                      par.mat = par.mat, par.clu = par.clu, Q.eff = Qeff.part[x])
    
    if((Qeff.part[x] == "par") || (Qeff.part[x] == "anc")){
      QTL <- QTL[, -1, drop = FALSE] }
    
  }
  
  QTL.el <- dim(QTL)[2] # number of QTL elements
  
  ### 2.1 cofactors
  
  cof.mat <- do.call(cbind, cof.list[which(cof.part[x, ])])
  
  # test if no cofactors
  
  if(is.null(cof.mat)){ cof.mat <- rep(0, length(mppData$geno.id)); cof.el <- 1
  
  } else { cof.el <- dim(cof.mat)[2] }
  
  
  # 2. model computation
  ######################
  
  ### 2.1 homogeneous residual variance error
  
  if(VCOV == "h.err"){
    
    model <- tryCatch(expr = lm(mppData$trait[, 1] ~ - 1 + cross.mat + cof.mat
                                + QTL), error = function(e) NULL)
    
    if (is.null(model)){
      
      results <- 0
      
    } else {
      
      
      if(!("QTL" %in% rownames(anova(model)))){ # QTL effect could not be
        # estimated probably due to
        # singularities.
        
        results <- 0
        
      } else {
        
        results <- -log10(anova(model)$Pr[which(rownames(anova(model))=="QTL")])
        
      }
      
    }
    
    ### 2.2 HRT REML or cross-specific variance residual terms
    
  } else if ((VCOV == "h.err.as") || (VCOV == "cr.err")){
    
    dataset <- data.frame(cof.mat = cof.mat, QTL = QTL,
                          cr.mat = factor(mppData$cross.ind,
                                          levels = unique(mppData$cross.ind)),
                          trait = mppData$trait[, 1])
    
    colnames(dataset) <- c(paste0("cof", 1:cof.el), paste0("Q", 1:QTL.el),
                           "cr.mat", "trait")
    
    formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
    formula.fix <- paste("trait ~ -1 + cr.mat + grp(cof)", formula.QTL)
    
    if(VCOV == "h.err.as"){ formula.R <- "~idv(units)"
    } else if (VCOV == "cr.err") {formula.R <- "~at(cr.mat):units"}
    
    model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                    rcov =  as.formula(formula.R),
                                    group = list(cof=1:cof.el),
                                    data=dataset, trace = FALSE,
                                    na.method.Y = "omit",
                                    na.method.X = "omit"),
                      error = function(e) NULL)
    
    ### 2.3 random pedigree + HVRT or + CSRT
    
  } else if ((VCOV == "pedigree") || (VCOV == "ped_cr.err")){
    
    dataset <- data.frame(cof.mat = cof.mat, QTL = QTL,
                          cr.mat = factor(mppData$cross.ind,
                                          levels = unique(mppData$cross.ind)),
                          trait = mppData$trait[, 1],
                          genotype = mppData$geno.id)
    
    colnames(dataset) <- c(paste0("cof", 1:cof.el), paste0("Q", 1:QTL.el),
                           "cr.mat", "trait", "genotype")
    
    
    formula.QTL <- paste("+",paste0("Q",1:QTL.el),collapse = " ")
    formula.fix <- paste("trait~1+grp(cof)",formula.QTL)
    
    if(VCOV == "pedigree"){ formula.R <- "~idv(units)"
    } else if (VCOV == "ped_cr.err") {formula.R <- "~at(cr.mat):units"}
    
    model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                    random = ~ ped(genotype),
                                    rcov = as.formula(formula.R),
                                    ginverse = list(genotype = ped.mat.inv),
                                    group = list(cof = 1:cof.el),
                                    data = dataset, trace = FALSE,
                                    na.method.Y = "omit",
                                    na.method.X = "omit"),
                      error = function(e) NULL)
    
  }
  
  # 3. organise the results for the mixed models similar for all models
  #####################################################################
  
  if(VCOV != "h.err"){
    
    if (is.null(model)){
      
      results <- 0
      
    } else {
      
      W.stat <- sum(wald(model)[3:(QTL.el+2), 3])
      
      if(W.stat == 0){
        
        results <- 0
        
      } else {
        
        df <- sum(wald(model)[3:(QTL.el+2), 3] != 0)
        
        pval <- pchisq(W.stat, df, lower.tail = FALSE)
        
        results <- -log10(pval)
        
      }
      
    }
    
  }
  
  return(results)
  
}