##########################
# MQE_eff_res_processing #
##########################

# function to process the computation of genetic effects from a MQE model

# x: vector indicating the QTL index

# mppData, mppData_bi : mppData objects

# Q.eff : character vector specifying the type of QTL effect

# results: result from a mixed or a linear model

# Q.pos: index of the QTL position

# Q.list: list of QTL incidence matrices

# par.clu: parent clustering object


MQE_eff_res_processing <- function(x, mppData, mppData_bi, Q.eff, results,
                                   Q.pos, Q.list, par.clu, VCOV){
  
  Qeff.i <-  Q.eff[x]
  res.i <-  results[[x]]
  Qpos.i <-  Q.pos[x]
  QTL.mat <-  Q.list[[x]]
  
  if(Qeff.i == "cr"){
    
    # projection: here no projection
    
    # add parents scores
    
    if(!is.null(mppData$geno.par)){
      
      par.sc_i <- mppData$geno.par[Qpos.i, 5:dim(mppData$geno.par)[2]]
      
      PA.all <- par.sc_i[mppData$par.per.cross[, 2]]
      PB.all <- par.sc_i[mppData$par.per.cross[, 3]]
      
      par.data <- data.frame(mppData$par.per.cross[, 2],
                             t(PA.all), mppData$par.per.cross[, 3],
                             t(PB.all), stringsAsFactors = FALSE)
      
      sign.eff <- sign(res.i[, 1])
      
      par.allele <- function(x, par.data, sign){
        
        if(is.na(sign[x])) {c(NA, NA)
        } else if(sign[x] < 0){par.data[x, 1:2]
        } else if (sign[x] > 0) {par.data[x, 3:4]}
        
      }
      
      par.sc_i <- lapply(X = 1:mppData$n.cr, FUN = par.allele,
                         sign = sign.eff, par.data = par.data)
      
      par.sc_i <- matrix(unlist(par.sc_i), ncol = 2, byrow = TRUE)
      
      
      results <- data.frame(res.i, par.sc_i, stringsAsFactors = FALSE)
      
      results[, 1] <- abs(results[, 1])
      
      colnames(results) <- c("Effect", "Std. Err.", "Test stat.", "p-value",
                             "sign.", "Add. parent", "All. parent")
      
      if(VCOV == "h.err"){colnames(results)[3] <- "t-test"} else {
        colnames(results)[3] <- "W-stat"}
      
      return(results)
      
    } else {
      
      return(res.i)
      
    }
    
    
  } else if (Qeff.i == "par"){
    
    # fill the reference parent value
    
    ref.mat <- matrix(0, nrow = mppData$n.par, ncol = 3)
    ref.mat <- cbind(ref.mat, rep(1, mppData$n.par))
    
    ref.names <- paste0("Q", x, mppData$parents)
    ref.mat[ref.names %in% rownames(res.i), ] <- as.matrix(res.i[, -5])
    rownames(ref.mat) <- ref.names
    
    # add sign stars back
    
    ref.mat <- data.frame(ref.mat, sapply(ref.mat[, 4], FUN = sign.star),
                          stringsAsFactors = FALSE)
    
    colnames(ref.mat) <- c("Effect", "Std. Err.", "Test stat.", "p-value",
                           "sign.")
    
    if(VCOV == "h.err"){colnames(ref.mat)[3] <- "t-test"} else {
      colnames(ref.mat)[3] <- "W-stat"}
    
    # addition of parents scores
    
    if(!is.null(mppData$geno.par)){
      
      par.sc_i <- mppData$geno.par[Qpos.i, 5:dim(mppData$geno.par)[2]]
      
      results <- data.frame(ref.mat, unlist(par.sc_i), stringsAsFactors = FALSE)
      colnames(results)[6] <- "All. parent"
      
      return(results)
      
    } else {
      
      return(ref.mat)
      
    }
    
    
    
  } else if (Qeff.i == "anc"){
    
    # fill the reference allele value
    
    n.anc.al <- dim(QTL.mat)[2] + 1 # number of ancestral allele
    
    if((n.anc.al == 2) && (VCOV == "h.err")){
      
      ref.mat <- matrix(0, nrow = n.anc.al, ncol = 3)
      ref.mat <- cbind(ref.mat, rep(1, n.anc.al))
      
      # reference alleles names
      
      allele <- sort(unique(par.clu[Qpos.i, ]))
      all.ind <- paste0("Q", x, "A.allele", allele)
      
      # QTL allele name
      
      Q.name <- paste0("Q", x, colnames(QTL.mat))
      
      ref.mat[all.ind %in% Q.name, ] <- as.matrix(res.i[, -5])
      rownames(ref.mat) <- all.ind
      
      
    } else {
      
      ref.mat <- matrix(0, nrow = n.anc.al, ncol = 3)
      ref.mat <- cbind(ref.mat, rep(1, n.anc.al))
      
      # reference names
      
      allele <- sort(unique(par.clu[Qpos.i, ]))
      all.ind <- paste0("Q", x, "A.allele", allele)
      
      ref.mat[all.ind %in% rownames(res.i), ] <- as.matrix(res.i[, -5])
      rownames(ref.mat) <- all.ind
      
    }
    
    
    # project into the parents (Attention NA values)
    
    A.allele <- factor(par.clu[Qpos.i, ])
    A <- model.matrix(~ A.allele - 1)
    ref.mat[is.na(ref.mat)] <- 9999 # replace NA before proejction
    ref.mat <- A %*% ref.mat
    ref.mat[ref.mat == 9999] <- NA
    
    # add sign stars back
    
    ref.mat <- data.frame(ref.mat, sapply(ref.mat[, 4], FUN = sign.star),
                          stringsAsFactors = FALSE)
    
    colnames(ref.mat) <- c("Effect", "Std. Err.", "Test stat.", "p-value",
                           "sign.")
    
    if(VCOV == "h.err"){colnames(ref.mat)[3] <- "t-test"} else {
      colnames(ref.mat)[3] <- "W-stat"}
    
    rownames(ref.mat) <- paste0("Q", x, mppData$parents)
    
    # add the parents values
    
    if(!is.null(mppData$geno.par)){
      
      par.sc_i <- mppData$geno.par[Qpos.i, 5:dim(mppData$geno.par)[2]]
      
      results <- data.frame(ref.mat, unlist(par.sc_i), stringsAsFactors = FALSE)
      colnames(results)[6] <- "All. parent"
      
      return(results)
      
    } else {
      
      return(ref.mat)
      
    }
    
    
  } else if (Qeff.i == "biall"){
    
    # projection: here no projection
    
    # add parents scores
    
    if(!is.null(mppData$geno.par)){
      
      ref.mat <- data.frame(matrix(0, nrow = mppData$n.par, ncol = 3),
                            rep(1, mppData$n.par), rep("", mppData$n.par),
                            stringsAsFactors = FALSE)
      
      par.sc_i <- mppData_bi$geno.par[Qpos.i, 5:dim(mppData$geno.par)[2]]
      
      ref.all <- mppData_bi$allele.ref[1, Qpos.i]
      het.sc <- mppData_bi$allele.ref[c(3, 4), Qpos.i]
      
      ind.na <- which(is.na(par.sc_i))
      ind.ref <- which(par.sc_i == ref.all)
      ind.het <- which(((par.sc_i == het.sc[1])|(par.sc_i == het.sc[2])))
      
      ref.mat[ind.ref, ] <- res.i
      ref.mat[ind.na, ] <- NA
      ref.mat[ind.het, 1] <- res.i[, 1]/2
      
      ref.mat <- data.frame(ref.mat, unlist(par.sc_i), stringsAsFactors = FALSE)
      
      rownames(ref.mat) <- paste0("Q", x, mppData_bi$parents)
      colnames(ref.mat) <- c("Effect", "Std. Err.", "Test stat.", "p-value",
                             "sign.", "All. parent")
      
      if(VCOV == "h.err"){colnames(ref.mat)[3] <- "t-test"} else {
        colnames(ref.mat)[3] <- "W-stat"}
      
      return(ref.mat)
      
    } else {
      
      return(res.i)
      
    }
    
  }
  
}