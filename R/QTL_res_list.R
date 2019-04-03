################
# QTL_res_list #
################

#' List of QTL results
#' 
#' Form a list of QTL results appending QTL effects results obtained during
#' different QTL detection procedure. The results are appended to
#' \code{res_file}.
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param MPP_out Output from \code{\link{mpp_proc}}.
#' 
#' @param trait \code{character} indicator to specify the  trait name.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 'cr', 'par', 'anc', 'biall' or 'MQE'.
#' 
#' @param VCOV \code{Character} expression defining the type of variance
#' covariance structure used.
#' 
#' @param res_file \code{data.frame} to store the QTL effects results. Default,
#' empty file.
#' 
#' @return 
#' 
#' The results of \code{MPP_out} are appended to \code{res_file}.
#' 
#' @author Vincent Garin
#' 
#' @examples
#' 
#' # not yet
#' 
#' @export
#'

QTL_res_list <- function(mppData, MPP_out, trait, Q.eff, VCOV,
                         res_file = c()){
  
  if(Q.eff == 'cr'){
    
    for(i in 1:MPP_out$n.QTL){
      
      cr_ind <- mppData$par.per.cross[, 1]
      cr_ind <- rep(cr_ind, each = 2)
      
      # Define the additive parents
      
      add_ind <- (MPP_out$QTL.effects[[1]][[i]]$Add.parent ==
                    mppData$par.per.cross[, 3]) * 1
      
      add_ind2 <- c()
      par_ind <- c()
      Qcr_eff_i <- c()
      
      for(j in 1:mppData$n.cr){
        
        Qcr_eff_i <- rbind(Qcr_eff_i, MPP_out$QTL.effects[[1]][[i]][j, ],
                           MPP_out$QTL.effects[[1]][[i]][j, ])
        
        if(add_ind[j] == 0 || is.na(add_ind[j])){
          
          add_ind2 <- c(add_ind2, c(1, -1))
          par_ind <- c(par_ind, mppData$par.per.cross[j, 2:3])
          
        } else {
          
          add_ind2 <- c(add_ind2, c(-1, 1))
          par_ind <- c(par_ind, rev(mppData$par.per.cross[j, 2:3]))
          
        }
        
      }
      
      Qcr_eff_i[, 1] <- Qcr_eff_i[, 1] * add_ind2
      
      Q_res_i <- data.frame(Proc = 'CIM', Trait = trait, Q.eff = Q.eff,
                            VCOV = VCOV, QTL_nb = i, MPP_out$QTL[i, -3],
                            MPP_out$QTL.CI[i, 4:8], cross = cr_ind,
                            parent = par_ind, Qcr_eff_i[, 1:5], Con.part = NA,
                            Par.all = NA, stringsAsFactors = FALSE,
                            row.names = NULL)
      
      res_file <- rbind(res_file, Q_res_i)
      
    }
    
    
    
  }
  
  if((Q.eff == 'par') || (Q.eff == 'anc')){
    
    for(i in 1:MPP_out$n.QTL){
      
      Q_res_i <- data.frame(Proc = 'CIM', Trait = trait, Q.eff = Q.eff,
                            VCOV = VCOV, QTL_nb = i, MPP_out$QTL[i, -3],
                            MPP_out$QTL.CI[i, 4:8], cross = NA,
                            parent = rownames(MPP_out$QTL.effects[[1]][[i]]),
                            MPP_out$QTL.effects[[1]][[i]],
                            stringsAsFactors = FALSE, row.names = NULL)
      
      res_file <- rbind(res_file, Q_res_i)
      
    }
    
  }
  
  if(Q.eff == 'biall'){
    
    for(i in 1:MPP_out$n.QTL){
      
      Q_res_i <- data.frame(Proc = 'CIM', Trait = trait, Q.eff = Q.eff,
                            VCOV = VCOV, QTL_nb = i, MPP_out$QTL[i, -3],
                            MPP_out$QTL.CI[i, 4:8], cross = NA,
                            parent = rownames(MPP_out$QTL.effects[[1]][[i]]),
                            MPP_out$QTL.effects[[1]][[i]][, 1:5],
                            Con.part = NA,
                            Par.all = MPP_out$QTL.effects[[1]][[i]][, 6],
                            stringsAsFactors = FALSE, row.names = NULL)
      
      res_file <- rbind(res_file, Q_res_i)
      
    }
    
  }
  
  if(Q.eff == 'MQE'){
    
    for(i in 1:MPP_out$n.QTL){
      
      if(MPP_out$QTL$QTL.eff[i] == 'cr'){
        
        cr_ind <- mppData$par.per.cross[, 1]
        cr_ind <- rep(cr_ind, each = 2)
        
        # Define the additive parents
        
        add_ind <- (MPP_out$QTL.effects[[i]]$Add.parent ==
                      mppData$par.per.cross[, 3]) * 1
        
        add_ind2 <- c()
        par_ind <- c()
        Qcr_eff_i <- c()
        
        for(j in 1:mppData$n.cr){
          
          Qcr_eff_i <- rbind(Qcr_eff_i, MPP_out$QTL.effects[[i]][j, ],
                             MPP_out$QTL.effects[[i]][j, ])
          
          if(add_ind[j] == 0 || is.na(add_ind[j])){
            
            add_ind2 <- c(add_ind2, c(1, -1))
            par_ind <- c(par_ind, mppData$par.per.cross[j, 2:3])
            
          } else {
            
            add_ind2 <- c(add_ind2, c(-1, 1))
            par_ind <- c(par_ind, rev(mppData$par.per.cross[j, 2:3]))
            
          }
          
        }
        
        Qcr_eff_i[, 1] <- Qcr_eff_i[, 1] * add_ind2
        
        range_info <- matrix(NA, ncol = 5)
        colnames(range_info) <- c('inf.mk', 'inf.lim[cM]', 'sup.mk',
                                  'sup.lim[cM]', 'range[cM]')
        
        Q_res_i <- data.frame(Proc = 'MQE', Trait = trait,
                              Q.eff = MPP_out$QTL$QTL.eff[i],
                              VCOV = VCOV, QTL_nb = i, MPP_out$QTL[i, -c(3, 5)],
                              log10pval = NA, range_info, cross = cr_ind,
                              parent = par_ind, Qcr_eff_i[, 1:5],
                              Con.part = NA, Par.all = NA,
                              stringsAsFactors = FALSE, row.names = NULL)
        
        res_file <- rbind(res_file, Q_res_i)
        
      }
      
      if((MPP_out$QTL$QTL.eff[i] == 'par')|| (MPP_out$QTL$QTL.eff[i] == 'anc')){
        
        range_info <- matrix(NA, ncol = 5)
        colnames(range_info) <- c('inf.mk', 'inf.lim[cM]', 'sup.mk',
                                  'sup.lim[cM]', 'range[cM]')
        
        Q_res_i <- data.frame(Proc = 'MQE', Trait = trait,
                              Q.eff = MPP_out$QTL$QTL.eff[i],
                              VCOV = VCOV, QTL_nb = i, MPP_out$QTL[i, -c(3, 5)],
                              log10pval = NA, range_info, cross = NA,
                              parent = rownames(MPP_out$QTL.effects[[i]]),
                              MPP_out$QTL.effects[[i]],
                              stringsAsFactors = FALSE, row.names = NULL)
        
        res_file <- rbind(res_file, Q_res_i)
        
      }
      
      
      if(MPP_out$QTL$QTL.eff[i] == 'biall'){
        
        range_info <- matrix(NA, ncol = 5)
        colnames(range_info) <- c('inf.mk', 'inf.lim[cM]', 'sup.mk',
                                  'sup.lim[cM]', 'range[cM]')
        
        Q_res_i <- data.frame(Proc = 'MQE', Trait = trait,
                              Q.eff = MPP_out$QTL$QTL.eff[i],
                              VCOV = VCOV, QTL_nb = i, MPP_out$QTL[i, -c(3, 5)],
                              log10pval = NA, range_info, cross = NA,
                              parent = rownames(MPP_out$QTL.effects[[i]]),
                              MPP_out$QTL.effects[[i]][, 1:5],
                              Con.part = NA,
                              Par.all = MPP_out$QTL.effects[[i]][, 6],
                              stringsAsFactors = FALSE, row.names = NULL)
        
        res_file <- rbind(res_file, Q_res_i)
        
      }
      
    }
    
  }
  
  return(res_file)
  
}