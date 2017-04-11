################
# QTLModelSIM #
################

# function to compute a single position QTL model

QTLModelSIM <- function(x, mppData, cross.mat, par.mat, Q.eff, par.clu, VCOV,
                        est.gen.eff){
  
  # 1. formation of the QTL incidence matrix
  ###########################################
  
  QTL <- IncMat_QTL(x = x, mppData = mppData, cross.mat = cross.mat,
                    par.mat = par.mat, par.clu = par.clu, Q.eff = Q.eff,
                    order.MAF = TRUE)
  
  QTL.el <- dim(QTL)[2] # number of QTL elements
  
  ref.name <- colnames(QTL)
  
  # 2. model computation
  ######################
  
  ### 2.1 homogeneous residual variance error
  
  if(VCOV == "h.err"){
    
    model <- tryCatch(expr = lm(mppData$trait[, 1] ~ - 1 + cross.mat + QTL),
                      error = function(e) NULL)
    
    if (is.null(model)){ 
      
      if(est.gen.eff) {
        
        if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
        
        } else { results <- c(0, rep(1, mppData$n.par)) }
        
      } else { results <- 0 }
      
    } else {
      
      if(!("QTL" %in% rownames(anova(model)))){ # QTL effect could not be
        # estimated probably due to
        # singularities.
        
        if(est.gen.eff) {
          
          if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
          
          } else { results <- c(0, rep(1, mppData$n.par)) }
          
        } else { results <- 0 }
        
      } else {
        
        if(est.gen.eff){
          
          gen.eff <- QTL_pval(mppData = mppData, model = model,
                              Q.eff = Q.eff, x = x, par.clu = par.clu)
          
          results <- c(-log10(anova(model)$Pr[2]), gen.eff)
          
        } else { results <- -log10(anova(model)$Pr[2]) }
        
      }
      
    }
    
    ### 2.2 HRT REML or cross-specific variance residual terms
    
  } else if ((VCOV == "h.err.as") || (VCOV == "cr.err")){
    
    dataset <- data.frame(QTL = QTL,
                          cr.mat = factor(mppData$cross.ind,
                                          levels = unique(mppData$cross.ind)),
                          trait = mppData$trait[, 1])
    colnames(dataset)[1:QTL.el] <- paste0("Q", 1:QTL.el)
    
    formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
    formula.fix <- paste("trait~-1+cr.mat", formula.QTL)
    
    if(VCOV == "h.err.as"){ formula.R <- "~idv(units)"
    } else if (VCOV == "cr.err") {formula.R <- "~at(cr.mat):units"} 
    
    model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                    rcov =  as.formula(formula.R),
                                    data = dataset, trace = FALSE,
                                    na.method.Y = "omit",
                                    na.method.X = "omit"),
                      error = function(e) NULL)
    
    
    
    ### 2.3 random pedigree + HVRT or + CSRT
    
  } else if ((VCOV == "pedigree") || (VCOV == "ped_cr.err")){
    
    dataset <- data.frame(QTL = QTL,
                          cr.mat = factor(mppData$cross.ind,
                                          levels = unique(mppData$cross.ind)),
                          trait = mppData$trait[, 1], genotype = mppData$geno.id)
    colnames(dataset)[1:QTL.el] <- paste0("Q", 1:QTL.el)
    
    
    formula.QTL <- paste("+", paste0("Q", 1:QTL.el), collapse = " ")
    formula.fix <- paste("trait~1", formula.QTL)
    
    if(VCOV == "pedigree"){ formula.R <- "~idv(units)"
    } else if (VCOV == "ped_cr.err") {formula.R <- "~at(cr.mat):units"}
    
    model <- tryCatch(expr = asreml(fixed = as.formula(formula.fix),
                                    random = ~ ped(genotype),
                                    rcov = as.formula(formula.R),
                                    ginverse = list(genotype = ped.mat.inv),
                                    data = dataset, trace = FALSE,
                                    na.method.Y = "omit",
                                    na.method.X = "omit"),
                      error = function(e) NULL)
    
  } 
  
  # 3. organise the results for the mixed models similar for all models
  #####################################################################
  
  if(VCOV != "h.err"){
    
    if (is.null(model)){ 
      
      if(est.gen.eff) {
        
        if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
        
        } else { results <- c(0, rep(1, mppData$n.par)) }
        
      } else { results <- 0 }
      
    } else {
      
      W.stat <- sum(wald(model)[2:(QTL.el+1), 3])
      
      if(W.stat == 0){
        
        if(est.gen.eff) {
          
          if(Q.eff == "cr"){ results <- c(0, rep(1, mppData$n.cr))
          
          } else { results <- c(0, rep(1, mppData$n.par)) }
          
        } else { results <- 0 }
        
      } else {
        
        df <- sum(wald(model)[2:(QTL.el+1), 1])
        
        pval <- pchisq(W.stat, df, lower.tail = FALSE)
        
        results <- -log10(pval)
        
        if(est.gen.eff){
          
          gen.eff  <- QTL_pval_mix(model = model, Q.eff = Q.eff, QTL.el = QTL.el,
                                   x = x, par.clu = par.clu,
                                   ref.name = ref.name, fct = "SIM")
          
          results  <- c(results, gen.eff)
          
        }
        
      }
      
    }
    
  }
  
  return(results)
  
}