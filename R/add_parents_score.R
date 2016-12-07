#####################
# add_parents_score #
#####################

# function to add the parents scores to the QTL effects results

add_parents_score <- function(results, mppData, Q.eff, QTL, par.ref, VCOV){
  
  n.QTL <- dim(QTL)[1]
  
  if (Q.eff == "cr"){
    
    Q.index <- paste0("Q", 1:n.QTL)
    Q.ind <- rep(Q.index, each = mppData$n.cr)
    
    parents.results <- c()
    
    for(i in 1:n.QTL){
      
      QTL.info <- results[Q.ind == Q.index[i], ]
      
      par.score <- mppData$geno.par[mppData$geno.par[, 1] == QTL[i, 1],
                                    5:dim(mppData$geno.par)[2]]
      
      PA.all <- par.score[mppData$par.per.cross[, 2]]
      PB.all <- par.score[mppData$par.per.cross[, 3]]
      
      par.data <- data.frame(mppData$par.per.cross[, 2],
                             t(PA.all), mppData$par.per.cross[, 3],
                             t(PB.all), stringsAsFactors = FALSE)
      
      sign.eff <- sign(QTL.info[, 1])
      
      par.allele <- function(x, par.data, sign){
        
        if(is.na(sign[x])) {c(NA, NA)
        } else if(sign[x] < 0){par.data[x, 1:2]
        } else if (sign[x] > 0) {par.data[x, 3:4]}
        
      }
      
      par.sc_i <- lapply(X = 1:mppData$n.cr, FUN = par.allele,
                         sign = sign.eff, par.data = par.data)
      
      par.sc_i <- matrix(unlist(par.sc_i), ncol = 2, byrow = TRUE)
      
      
      parents.results <- rbind(parents.results, par.sc_i)
      
    }
    
    results[, 1] <- abs(results[, 1])
    results <- data.frame(results, parents.results,
                          stringsAsFactors = FALSE)
    
    colnames(results) <- c("Effect", "Std. Err.", "Test stat.", "p-value",
                           "sign.", "Add. parent", "All. parent")
    
    if(VCOV == "h.err"){colnames(results)[3] <- "t-test"} else {
      colnames(results)[3] <- "W-stat"}
    
  } else if ((Q.eff == "par") || (Q.eff == "anc")){
    
    
    parents.results <- c()
    
    for(i in 1:n.QTL){
      
      par.sc_i <- mppData$geno.par[mppData$geno.par[, 1] == QTL[i, 1],
                                   5:dim(mppData$geno.par)[2]]
      
      parents.results <- c(parents.results, unlist(par.sc_i))
      
    }
    
    results <- data.frame(results, parents.results,
                          stringsAsFactors = FALSE)
    
    colnames(results) <- c("Effect", "Std. Err.", "Test stat.", "p-value",
                           "sign.", "All. parent")
    
    if(VCOV == "h.err"){colnames(results)[3] <- "t-test"} else {
      colnames(results)[3] <- "W-stat"}
    
  } else if (Q.eff == "biall"){
    
    # find the reference allele at the QTL positions
    
    parents.results <- c()
    
    for(i in 1:n.QTL){
      
      ref.mat <- data.frame(matrix(0, nrow = mppData$n.par, ncol = 3),
                            rep(1, mppData$n.par), rep("", mppData$n.par),
                            stringsAsFactors = FALSE)
      
      index <- which(mppData$geno.par[, 1] == QTL[i, 1])
      par_sc_Qi <- mppData$geno.par[index, 5:dim(mppData$geno.par)[2]]
      par_sc_Qi <- unlist(par_sc_Qi)
      
      
      ref.all <- c(mppData$allele.ref[1, index, drop = FALSE])
      
      het.sc <- mppData$allele.ref[c(3, 4), index]
      
      ind.na <- which(is.na(par_sc_Qi))
      ind.ref <- which(par_sc_Qi == ref.all)
      ind.het <- which(((par_sc_Qi == het.sc[1])|(par_sc_Qi == het.sc[2])))
      
      
      ref.mat[ind.ref, ] <- results[i, ]
      ref.mat[ind.na, ] <- NA
      ref.mat[ind.het, 1] <- results[i, 1]/2
      
      ref.mat <- data.frame(ref.mat, par_sc_Qi, stringsAsFactors = FALSE)
      
      rownames(ref.mat) <- paste0("Q",i, mppData$parents)
      colnames(ref.mat) <- c("Effect", "Std. Err.", "Test stat.", "p-value",
                             "sign.", "All. parent")
      
      if(VCOV == "h.err"){colnames(ref.mat)[3] <- "t-test"} else {
        colnames(ref.mat)[3] <- "W-stat"}
      
      parents.results <- rbind(parents.results, ref.mat)
      
    }
    
    results <- parents.results
    
  }
  
  return(results)
}