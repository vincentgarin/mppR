#######################
# Qeff_res_processing #
#######################

# Processing of the QTL results

# arguments

# model : multi-QTL model

# mppData : data object

# Q.list : list of QTLs

# cross.mat : cross intercept incidence matrix

# VCOV : 

# Q.eff

# QTL: list of QTLs

# par.clu: parental clustering object

# par.ref: parent reference indicator

# const: type of constraint

Qeff_res_processing <- function(model, mppData, cross.mat, Q.list, QTL,
                                Q.eff, par.clu, VCOV, par.ref, const){
  
  n.QTL <- length(Q.list)
  
  # deal with the two case where 
  
  if(const == "sum.0"){
    
    results <- summary(model)$coefficients
    index <- (substr(rownames(results), 1, 1) == "Q")
    results <- subset(x = results, subset = index, drop = FALSE)
    
    if(Q.eff == "par"){
      
      # fill the reference parent value
      
      ref.mat <- matrix(NA, nrow = (mppData$n.par* n.QTL), ncol = 4)
      
      # fill the potential singular values
      
      parents.ind <- rep(mppData$parents, times = n.QTL)
      ref.names <- paste0(rep(paste0("Q", 1:n.QTL), each = mppData$n.par),
                          parents.ind)
      ref.mat[ref.names %in% rownames(results), ] <- results
      rownames(ref.mat) <- ref.names
      
      return(ref.mat)
      
    } else if (Q.eff == "anc"){
      
      
      n.anc.al <- unlist(lapply(X = Q.list, FUN = function(x) dim(x)[2]))
      Q.names <- names(n.anc.al)
      QTL.ind <- rep(Q.names, time = n.anc.al)
      n.anc.al <- sum(n.anc.al)
      
      ref.mat <- matrix(NA, nrow = n.anc.al, ncol = 4)
      
      anc.all.name <- function(x, Q.list) {paste0("Q", x, colnames(Q.list[[x]]))}
      
      all.ind <- unlist(lapply(X = 1:n.QTL, FUN = anc.all.name, Q.list = Q.list))
      
      
      ref.mat[(all.ind %in% rownames(results)), ] <- results
      rownames(ref.mat) <- all.ind
      
      # project into the parents
      
      Qeff.mat <- c()
      
      for(i in 1:n.QTL){
        
        Q.mat <- ref.mat[QTL.ind == names(Q.list)[i], ]
        A.allele <- factor(par.clu[which(rownames(par.clu) == QTL[i, 1]), ])
        A <- model.matrix(~ A.allele - 1)
        Q.mat[is.na(Q.mat)] <- 9999 # replace NA before proejction
        Q.mat <- A %*% Q.mat
        Q.mat[Q.mat == 9999] <- NA
        Qeff.mat <- rbind(Qeff.mat, Q.mat)
        
      }
      
      colnames(Qeff.mat) <- colnames(ref.mat)
      parents.ind <- rep(mppData$parents, n.QTL)
      ref.names <- paste0(rep(paste0("Q", 1:n.QTL), each = mppData$n.par),
                          parents.ind)
      rownames(Qeff.mat) <- ref.names
      
      return(Qeff.mat)
      
    }
    
  } else {
    
    if(VCOV == "h.err"){
      
      results <- summary(model)$coefficients
      index <- (substr(rownames(results), 1, 1) == "Q")
      results <- subset(x = results, subset = index, drop = FALSE)
      
      # control for singular values and fill missing values
      
      if (Q.eff == "cr"){
        
        ref.mat <- matrix(NA, nrow = (mppData$n.cr * n.QTL), ncol = 4)
        ref.names <- paste0(rep(paste0("Q", 1:n.QTL), each = mppData$n.cr),
                            rep(colnames(cross.mat), times = n.QTL))
        ref.mat[ref.names %in% rownames(results), ] <- results
        rownames(ref.mat) <- ref.names
        
        return(ref.mat)
        
      } else if (Q.eff == "par"){
        
        # fill the reference parent value
        
        ref.mat <- matrix(NA, nrow = (mppData$n.par* n.QTL), ncol = 4)
        parents.ind <- rep(mppData$parents, times = n.QTL)
        ref.mat[parents.ind %in% par.ref, ] <- matrix(rep(c(0, 0, 0, 1), n.QTL),
                                                      nrow = n.QTL, byrow = TRUE)
        
        # fill the potential singular values
        
        ref.names <- paste0(rep(paste0("Q", 1:n.QTL), each = mppData$n.par),
                            parents.ind)
        ref.mat[ref.names %in% rownames(results), ] <- results
        rownames(ref.mat) <- ref.names
        
        return(ref.mat)
        
      } else if (Q.eff == "anc") {
        
        # fill the reference ancestral allele value
        
        n.anc.al <- unlist(lapply(X = Q.list, FUN = function(x) dim(x)[2]))
        
        # change the name of the individual ancestral class
        
        if(sum(n.anc.al == 1) > 0){
          
          QTL.name <- names(n.anc.al[n.anc.al == 1])
          QTL.nb <- which(n.anc.al == 1)
          new.names <- lapply(X = QTL.nb,
                              FUN = function(x, Q.list) colnames(Q.list[[x]]),
                              Q.list = Q.list)
          new.names <- paste0(QTL.name, unlist(new.names))
          
          row.names <- rownames(results)
          row.names[row.names %in% QTL.name] <- new.names
          rownames(results) <- row.names
          
        }
        
        Q.names <- names(n.anc.al)
        QTL.ind <- rep(Q.names, time = (n.anc.al + 1))
        n.anc.al <- sum((n.anc.al + 1))
        
        ref.mat <- matrix(NA, nrow = n.anc.al, ncol = 4)
        
        anc.allele <- function(x, par.clu, QTL, par.ref){
          
          clu.info <- par.clu[which(rownames(par.clu) == QTL[x, 1]), ]
          clu.to.remove <- clu.info[names(clu.info) == par.ref]
          allele <- sort(unique(clu.info))
          
          ref.all <- paste0("Q", x, "A.allele", clu.to.remove)
          all.ind <- paste0("Q", x, "A.allele", allele)
          
          list(ref.all = ref.all, all.ind = all.ind)
          
        }
        
        anc.ind <- lapply(X = 1:n.QTL, FUN = anc.allele, par.clu = par.clu,
                          QTL = QTL, par.ref = par.ref)
        
        all.ind <- unlist(lapply(X = anc.ind, FUN = function(x) x$all.ind))
        ref.all <- unlist(lapply(X = anc.ind, FUN = function(x) x$ref.all))
        
        
        ref.mat[(all.ind %in% ref.all), ] <- matrix(rep(c(0, 0, 0, 1), n.QTL),
                                                    nrow = n.QTL, byrow = TRUE)
        
        # fill the potential singular values
        
        ref.mat[(all.ind %in% rownames(results)), ] <- results
        rownames(ref.mat) <- all.ind
        
        # project into the parents
        
        Qeff.mat <- c()
        
        for(i in 1:n.QTL){
          
          Q.mat <- ref.mat[QTL.ind == names(Q.list)[i], ]
          A.allele <- factor(par.clu[which(rownames(par.clu) == QTL[i, 1]), ])
          A <- model.matrix(~ A.allele - 1) 
          Q.mat[is.na(Q.mat)] <- 9999 # replace NA before proejction
          Q.mat <- A %*% Q.mat
          Q.mat[Q.mat == 9999] <- NA
          Qeff.mat <- rbind(Qeff.mat, Q.mat)
          
        }
        
        colnames(Qeff.mat) <- colnames(ref.mat)
        parents.ind <- rep(mppData$parents, n.QTL)
        ref.names <- paste0(rep(paste0("Q", 1:n.QTL), each = mppData$n.par),
                            parents.ind)
        rownames(Qeff.mat) <- ref.names
        
        return(Qeff.mat)
        
      } else if (Q.eff == "biall"){
        
        ref.mat <- matrix(NA, nrow = n.QTL, ncol = 4)
        ref.names <- paste0("Q", 1:n.QTL)
        
        ref.mat[ref.names %in% rownames(results), ] <- results
        rownames(ref.mat) <- ref.names
        
        return(ref.mat)
        
      }
      
      
      
    } else {
      
      # form the table with the results
      
      index <- substr(names(rev(model$coefficients$fixed)), 1, 1) == "Q"
      
      w.table <- wald(model)
      w.stat <- w.table[substr(rownames(w.table), 1, 1) == "Q", c(3, 4)]
      
      results <- data.frame(rev(model$coefficients$fixed)[index],
                            rev(sqrt(model$vcoeff$fixed))[index],
                            w.stat)
      
      # put as NA the singular values
      
      results[results[, 1] == 0, ] <- NA
      
      if((Q.eff == "cr") || (Q.eff == "biall")){
        
        return(results)
        
      } else if (Q.eff == "par"){
        
        # put the reference parent
        
        ref.mat <- matrix(0, nrow = (mppData$n.par* n.QTL), ncol = 4)
        ref.mat[, 4] <- 1
        
        Q.names <- rep(paste0("Q", 1:n.QTL), each = mppData$n.par)
        parents.ind <- rep(mppData$parents, times = n.QTL)
        parents.ind <- paste0(Q.names, parents.ind)
        
        ref.mat[parents.ind %in% rownames(results) , ] <- as.matrix(results)
        
        rownames(ref.mat) <- parents.ind
        
        return(ref.mat)
        
      } else if (Q.eff == "anc") {
        
        # fill the reference ancestral allele value
        
        n.anc.al <- unlist(lapply(X = Q.list, FUN = function(x) dim(x)[2]))
        
        Q.names <- names(n.anc.al)
        
        QTL.ind <- unlist(mapply(FUN = function(x, names) rep(names, (x+1)),
                                 x = n.anc.al, names = Q.names), use.names = FALSE)
        
        n.anc.al <- sum((n.anc.al + 1))
        
        ref.mat <- matrix(NA, nrow = n.anc.al, ncol = 4)
        
        anc.allele <- function(x, par.clu, QTL, par.ref){
          
          clu.info <- par.clu[which(rownames(par.clu) == QTL[x, 1]), ]
          clu.to.remove <- clu.info[names(clu.info) == par.ref]
          allele <- sort(unique(clu.info))
          
          ref.all <- paste0("Q", x, "A.allele", clu.to.remove)
          all.ind <- paste0("Q", x, "A.allele", allele)
          
          list(ref.all = ref.all, all.ind = all.ind)
          
        }
        
        anc.ind <- lapply(X = 1:n.QTL, FUN = anc.allele, par.clu = par.clu,
                          QTL = QTL, par.ref = par.ref)
        
        all.ind <- unlist(lapply(X = anc.ind, FUN = function(x) x$all.ind))
        ref.all <- unlist(lapply(X = anc.ind, FUN = function(x) x$ref.all))
        
        
        ref.mat[(all.ind %in% ref.all), ] <- matrix(rep(c(0, 0, 0, 1),
                                                        n.QTL),
                                                    nrow = n.QTL, byrow = TRUE)
        
        # fill the potential singular values
        
        ref.mat[(all.ind %in% rownames(results)), ] <- as.matrix(results)
        rownames(ref.mat) <- all.ind
        
        # project into the parents
        
        Qeff.mat <- c()
        
        for(i in 1:n.QTL){
          
          Q.mat <- ref.mat[QTL.ind == names(Q.list)[i], ]
          A.allele <- factor(par.clu[which(rownames(par.clu) == QTL[i, 1]), ])
          A <- model.matrix(~ A.allele - 1) 
          Q.mat[is.na(Q.mat)] <- 9999 # replace NA before proejction
          Q.mat <- A %*% Q.mat
          Q.mat[Q.mat == 9999] <- NA
          Qeff.mat <- rbind(Qeff.mat, Q.mat)
          
        }
        
        parents.ind <- rep(mppData$parents, n.QTL)
        ref.names <- paste0(rep(paste0("Q", 1:n.QTL), each = mppData$n.par),
                            parents.ind)
        rownames(Qeff.mat) <- ref.names
        
        return(Qeff.mat)
        
      } 
      
      
    }
    
  }

}