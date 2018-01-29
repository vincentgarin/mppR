###############
# MQE_forward #
###############

#' Forward regression with different type of QTL effects
#' 
#' Determines a multi-QTL effect (MQE) model using a forward regression.
#' 
#' The included QTL position can have different type of QTL effects at different
#' loci. At each step (new added position), the function compute QTL profiles
#' using one type of QTL effect specified in \code{Q.eff} for the tested
#' position. Let us assume that the user want to allow the QTL to have a
#' parental, ancestral or bi-allelic effect. Therefore, the function will
#' calculate three profiles. Then it will select in each profile the most
#' significant position if a position is above the \code{threshold} value.
#' The position and its QTL effect that increase the most the global adjusted R
#' squared (\code{MQE_R2}) will be selected as the new QTL position and included
#' in the list of cofactors. The function continues to iterate until no position
#' is significant anymore.
#' 
#' \strong{WARNING!(1)} The computation of \code{MQE_forward()} function using
#' mixed models (all models with \code{VCOV} different than \code{"h.err"})
#' is technically possible but can be irrealistic in practice due to a reduced
#' computer power. Since a mixed model is computed at each single position it
#' can take a lot of time. From our estimation it can take between 20 to 50
#' times more time than for the linear model (HRT). If the number of detected
#' QTL is supposed to be small (until 5) it could still be feasible.
#' 
#' \strong{WARNING!(2)} The computation of random pedigree models
#' (\code{VCOV = "pedigree" and "ped_cr.err"}) can sometimes fail. This could be
#' due to singularities due to a strong correlation between the QTL term(s) and 
#' the polygenic term. This situation can appear in the parental model.
#' the error can also sometimes come from the \code{asreml()} function. From
#' our experience, in that case, trying to re-run the function one or two times
#' allow to obtain a result.
#'
#' @param mppData An IBD object of class \code{mppData}
#' See \code{\link{mppData_form}} for details. Default = NULL.
#'
#' @param mppData_bi Required IBS object of class \code{mppData} if the user
#' wants to allow QTLs with a bi-allelic effect. \strong{The list of marker must
#' be strictly the same as the one of \code{mppData}.} Default = NULL.
#' 
#' @param Q.eff \code{Character} vector of possible QTL effects the user want to
#' test. Elements of Q.eff can be "cr", "par", "anc" or "biall". For details
#' look at \code{\link{mpp_SIM}}.
#' 
#' @param par.clu Required argument if the user wants to allow QTLs with an
#' ancestral effect. \code{Interger matrix} representing the results of a parents genotypes
#' clustering. The columns represent the parental lines and the rows
#' the different markers or in between positions. \strong{The columns names must
#' be the same as the parents list of the mppData object. The rownames must be
#' the same as the map marker list of the mppData object.} At a particular
#' position, parents with the same value are assumed to inherit from the same
#' ancestor. for more details, see \code{\link{USNAM_parClu}} and
#' \code{\link{parent_cluster}}. Default = NULL.
#'
#' @param VCOV \code{Character} expression defining the type of variance
#' covariance structure used: 1) "h.err" for an homogeneous variance residual term
#' (HRT) linear model; 2) "h.err.as" for a HRT model fitted by REML using
#' \code{ASReml-R}; 3) "cr.err" for a cross-specific variance residual terms
#' (CSRT) model; 4) "pedigree" for a random pedigree term and HRT model;
#' and 5) "ped_cr.err" for random pedigree and CSRT model.
#' For more details see \code{\link{mpp_SIM}}. Default = "h.err".
#' 
#' @param threshold \code{Numeric} value representing the -log10(p-value) threshold
#' above which a position can be considered as significant. Default = 3.
#' 
#' @param window \code{Numeric} distance (cM) on the left and the right of a
#' cofactor position where it is not included in the model. Default = 20.
#' 
#' @param parallel \code{Logical} value specifying if the function should be
#' executed in parallel on multiple cores. To run function in parallel user must
#' provide cluster in the \code{cluster} argument. \strong{Parallelization is
#' only available for HRT (linear) models \code{VCOV = "h.err"}}.
#' Default = FALSE.
#'
#' @param cluster Cluster object obtained with the function \code{makeCluster()}
#' from the parallel package. Default = NULL.
#' 
#' @param verbose \code{Logical} value indicating if the steps of the
#' MQE_forward should be printed. It will not affect the printing of the other
#' functions called by \code{mpp_proc()}, especially the printing of
#' \code{asreml()}. Default = TRUE.
#' 
#' 
#' @return
#' 
#' \item{QTL.list }{\code{Data.frame} with six columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; 5) -log10(p-values); and 6) type of
#' QTL incidence matrix of the selected positions}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mpp_perm}},
#' \code{\link{mpp_SIM}},
#' \code{\link{MQE_R2}}, \code{\link{parent_cluster}},
#' \code{\link{USNAM_parClu}}
#'
#' @examples
#'
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' data(USNAM_parClu)
#' 
#' mppData <- USNAM_mppData
#' mppData_bi <- USNAM_mppData_bi
#' par.clu <- USNAM_parClu
#' 
#' 
#' # equalize the list of markers of the two mppData objects and par.clu
#' 
#' com.mk.list <- intersect(mppData$map$mk.names, mppData_bi$map$mk.names)
#' 
#' mppData <- mppData_subset(mppData = mppData, mk.list = com.mk.list)
#' mppData_bi <- mppData_subset(mppData = mppData_bi, mk.list = com.mk.list)
#' par.clu <- par.clu[rownames(par.clu) %in% com.mk.list, ]
#' 
#' 
#' QTL <- MQE_forward(mppData = mppData, mppData_bi = mppData_bi,
#'                    Q.eff = c("par", "anc", "biall"), par.clu = par.clu)
#'                    
#' R2 <- MQE_R2(mppData = mppData, QTL = QTL[, 1], Q.eff = QTL[, 5],
#' par.clu = par.clu)
#'
#' @export
#'


MQE_forward <- function(mppData = NULL, mppData_bi = NULL, Q.eff, par.clu = NULL,
                        VCOV = "h.err", threshold = 3, window = 20,
                        parallel = FALSE, cluster = NULL,
                        verbose = TRUE) {
  
  # 1. test data format
  #####################
  
  check.MQE(mppData = mppData, mppData_bi = mppData_bi, Q.eff = Q.eff,
            VCOV = VCOV, par.clu = par.clu, parallel = parallel,
            cluster = cluster, fct = "forward")
  
  # Sub-function for R squared computation
  
  R2.sg <-  function(mppData, mppData_bi, QTL, Q.eff, par.clu){
    
    MQE_R2(mppData = mppData, mppData_bi = mppData_bi, QTL = QTL,
           Q.eff = Q.eff, par.clu = par.clu, glb.only = TRUE)[[2]]
    
  }
  
  # form eventual pedigree term
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  # check parental clustering object
  
  if ("anc" %in% Q.eff) {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  # initialize list of selected QTL and list of type of effects.
  
  QTL.list <- c()
  QTL.eff <- c()
  pos.ind <- 1
  
  # 2. step 1: SIM
  ################
  
  # vector to store the candidate positions
  
  max.pval <- c()
  cand.pos <- c()
  
  
  for (i in 1:length(Q.eff)) {
    
    Q.eff_i <- Q.eff[i]
    
    if (Q.eff_i == "biall") {
      
      SIM <- mpp_SIM(mppData = mppData_bi, Q.eff = Q.eff_i,
                     par.clu = par.clu, VCOV = VCOV, parallel = parallel,
                     cluster = cluster)
      
    } else {
      
      SIM <- mpp_SIM(mppData = mppData, Q.eff = Q.eff_i,
                     par.clu = par.clu, VCOV = VCOV, parallel = parallel,
                     cluster = cluster)
      
    }
    
    # store the candidate position
    
    max.pval_i <- max(SIM$log10pval, na.rm = TRUE)
    
    max.pval <- c(max.pval, max.pval_i)
    
    if (max.pval_i > threshold) {
      
      cand.pos_i <- data.frame(mppData$map[which.max(SIM$log10pval), 1], Q.eff_i,
                               stringsAsFactors = FALSE)
      
      cand.pos <- rbind.data.frame(cand.pos, cand.pos_i)
      
    }
    
  }
  
  
  if(!is.null(cand.pos)) {colnames(cand.pos) <- c("mk.names", "Q.eff")}
  
  if (!is.null(cand.pos)) {
    
    if(verbose){
      
      cat("\n")
      cat(paste("position", pos.ind))
      cat("\n")
      
    }
    
    R2.pos1 <- mapply(FUN = R2.sg, QTL = cand.pos[, 1], Q.eff = cand.pos[, 2],
                      MoreArgs = list(mppData = mppData, mppData_bi = mppData_bi,
                                      par.clu = par.clu))
    
    # test which position has the highest R squared
    
    if ((length(R2.pos1) != 1) && (length(unique(R2.pos1)) == 1)) {
      
      sel.pos_i <- cand.pos[which.max(max.pval), ]
      
      # two pos with max values, take the first one.
      
      if (dim(sel.pos_i)[1] > 1) { sel.pos_i <- sel.pos_i[1, ] }
      
    } else {
      
      sel.pos_i <- cand.pos[which.max(R2.pos1), ]
      
      if (dim(sel.pos_i)[1] > 1) {sel.pos_i <- sel.pos_i[1, ] }
      
    }
    
    QTL.list <- c(QTL.list, unlist(sel.pos_i[1]))
    QTL.eff <- c(QTL.eff, unlist(sel.pos_i[2]))
    pos.ind <- pos.ind + 1
    
    
  } else {
    
    
    warning(paste("No position is above the threshold, the stepwise procedure",
               "could not select any QTL position."))
    
    return(NULL)
    
  }
  
  # 3. Step 2 : CIM with cofactors (QTL already selected)
  #######################################################
  
  rest.sign.pos <- TRUE
  
  while(rest.sign.pos) {
    
    
    max.pval <- c()
    cand.pos <- c()
    
    for (i in 1:length(Q.eff)) {
      
      Q.eff_i <- Q.eff[i]
      
      log10.pval <- MQE_CIM(mppData = mppData, mppData_bi = mppData_bi,
                            Q.eff = Q.eff_i, par.clu = par.clu, VCOV = VCOV,
                            cofactors = QTL.list, cof.Qeff = QTL.eff,
                            parallel = parallel, cluster = cluster)
      
      
      # select the position with the highest p-value that is not in an
      # already selected QTL region.
      
      # make an index of positions that are not in a QTL region
      
      pos.QTL <- mppData$map[mppData$map[, 1] %in% QTL.list, c(2,4)]
      
      test.cof <- function(x, map, window) {
        
        t1 <- map$chr == as.numeric(x[1])
        t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
        (t1 & t2)
        
      }
      
      cof.part <- apply(X = pos.QTL, MARGIN = 1, FUN = test.cof,
                        map = mppData$map, window = window)
      
      log10.pval <- log10.pval[rowSums(cof.part) == 0, ]
      
      max.pval_i <- max(log10.pval$log10pval, na.rm = TRUE)
      
      
      if (max.pval_i > threshold) {
        
        max.pval <- c(max.pval, max.pval_i)
        
        cand.pos_i <- data.frame(log10.pval[which.max(log10.pval$log10pval), 1],
                                 Q.eff_i, stringsAsFactors = FALSE)
        
        cand.pos <- rbind.data.frame(cand.pos, cand.pos_i)
        
      }
      
    }
    
    if(!is.null(cand.pos)) {colnames(cand.pos) <- c("mk.names", "Q.eff")}
    
    # test if some position has been added. If yes, determine which one
    # increase the R squared
    
    if (!is.null(cand.pos)) {
      
      if(verbose){
        
        cat("\n")
        cat(paste("position", pos.ind))
        cat("\n")
        
      }
      
      # evaluate which position give the best extra adjusted R squared
      
      QTL.list_i <- lapply(X = cand.pos[, 1],
                           FUN = function(x, Q.l) c(Q.l, x), Q.l = QTL.list)
      
      QTL.eff_i <- lapply(X = cand.pos[, 2],
                          FUN = function(x, Q.e) c(Q.e, x), Q.e = QTL.eff)
      
      # R2 ith position
      
      R2.posi <- mapply(FUN = R2.sg, QTL = QTL.list_i, Q.eff = QTL.eff_i,
                        MoreArgs = list(mppData = mppData, mppData_bi = mppData_bi,
                                        par.clu = par.clu))
      
      # test which position has the highest R squared
      
      if ((length(R2.posi) != 1) & (length(unique(R2.posi)) == 1)) {
        
        sel.pos_i <- cand.pos[which.max(max.pval), ]
        
        if (dim(sel.pos_i)[1] > 1) { sel.pos_i <- sel.pos_i[1, ] }
        
      } else {
        
        sel.pos_i <- cand.pos[which.max(R2.posi), ]
        
        if (dim(sel.pos_i)[1] > 1) {sel.pos_i <- sel.pos_i[1, ] }
        
      }
      
      QTL.list <- c(QTL.list, unlist(sel.pos_i[1]))
      QTL.eff <- c(QTL.eff, unlist(sel.pos_i[2]))
      pos.ind <- pos.ind + 1
      
      
    } else { rest.sign.pos <- FALSE }
    
    
  }  # end while loop
  
  
  # return the final list of cofactors (QTL)
  
  
  results <- mppData$map[mppData$map[, 1] %in% QTL.list, ]
  names(QTL.eff) <- QTL.list
  QTL.eff2 <- QTL.eff[results[, 1]]
  results <- data.frame(results, QTL.eff2, stringsAsFactors = FALSE)
  
  return(results)
  
  
}