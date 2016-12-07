################
# MQE_BackElim #
################

#' Backward elimination on multi-QTL effect candidates
#' 
#' Performs a backward elimination using a list of given QTLs positions. These
#' position can have different type of QTL effects (cross-specific, parental,
#' ancestral or bi-allelic). The positions with a p-value above the significance
#' level \code{alpha}, are successively removed.
#' 
#' The function starts with all QTL positions in the model and test the inclusion
#' of each position as the last in the model. If all position p-value are below
#' \code{alpha} the procedure stop. If not the position with the highest p-value
#' is remove and the procedure continue until there are no more position with 
#' a p-value below \code{alpha}.
#' 
#' \strong{WARNING!} The estimation of the random pedigree models
#' (\code{VCOV = "pedigree" and "ped_cr.err"}) can be unstable. Sometimes the
#' \code{asreml()} function fails to produce a results and returns the following
#' message: \strong{\code{GIV matrix not positive definite: Singular pivots}}.
#' So far we were not able to identify the reason of this problem and to
#' reproduce this error because it seems to happen randomly. From our
#' experience, trying to re-run the function one or two times should allow
#' to obtain a result.
#' 
#' @param mppData An IBD object of class \code{mppData}
#' See \code{\link{mppData_form}} for details. Default = NULL.
#'
#' @param mppData_bi Required IBS object of class \code{mppData} if the user
#' wants to allow QTLs with a bi-allelic effect. \strong{The list of marker must
#' be strictly the same as the one of \code{mppData}.} Default = NULL.
#' 
#' @param QTL Vector of \code{character} markers or inbetween marker positions
#' names. Default = NULL.
#' 
#' @param Q.eff \code{Character} vector indicating for each QTL position the
#' type of QTL effect among: "cr", "par", "anc" and "biall". For details look at
#' \code{\link{mpp_SIM}}.
#'
#' @param par.clu Required argument if the user wants to allow QTLs with an
#' ancestral effect. \code{interger matrix} representing the results of a parents
#' genotypes
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
#' @param alpha \code{Numeric} value indicating the level of significance for
#' the backward elimination. Default = 0.05.
#' 
#' @return Return:
#' 
#' \item{QTL }{Object of class \code{QTLlist}. It is a \code{data.frame}
#' with four columns containing information on the QTL candidates : marker or
#' added position names, chromosome indicator, position in centi-Morgan,
#' -log10(p-val) and the type of QTL effect.}
#' 
#' \item{QTL }{\code{Data.frame} of class \code{QTLlist} with six columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) Interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; 5) -log10(p-values), and 6) type of QTL effect.}
#'  
#' @author Vincent Garin
#' 
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
#' # Equalize the list of markers of the two mppData objects and par.clu
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
#' QTL <-  MQE_BackElim(mppData = mppData, mppData_bi = mppData_bi,
#'                      QTL = QTL[, 1], Q.eff = QTL[, 5], par.clu = par.clu)
#' 
#' @export
#' 


MQE_BackElim <- function (mppData = NULL, mppData_bi = NULL, QTL = NULL, Q.eff,
                          par.clu = NULL, VCOV = "h.err", alpha = 0.05) {
  
  # 1. Check data format
  ######################
  
  check.MQE(mppData = mppData, mppData_bi =  mppData_bi, Q.eff = Q.eff,
            VCOV = VCOV, par.clu = par.clu, QTL = QTL, fct = "QTLeffects")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 inverse of the pedigree matrix
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.2 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.3 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.4 modify the par.clu object order parents columns and replace monomorphic
  
  if ("anc" %in% Q.eff) {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.5 Formation of the list of QTL
  
  # order list of QTL positions
  
  Q.pos <- vapply(X = QTL,
                  FUN = function(x, mppData) which(mppData$map[,1] == x),
                  FUN.VALUE = numeric(1), mppData = mppData)
  
  Q.ord <- data.frame(QTL, Q.eff, Q.pos, stringsAsFactors = FALSE)
  
  Q.ord <- Q.ord[order(Q.pos),]
  
  QTL <- Q.ord[, 1]; Q.eff <- Q.ord[, 2]
  
  Q.pos <- which(mppData$map[, 1] %in% QTL)
  
  
  # function to produce different type of QTL incidence matricdes
  
  IncMat_QTL_MQE <- function(x, mppData, mppData_bi, Q.eff, par.clu,
                             cross.mat, par.mat){
    
    if(Q.eff == "biall") {
      
      IncMat_QTL(x = x, mppData = mppData_bi, Q.eff = Q.eff, par.clu = par.clu,
                 cross.mat = cross.mat, par.mat = par.mat)
      
    } else {
      
      Q <- IncMat_QTL(x = x, mppData = mppData, Q.eff = Q.eff,
                      par.clu = par.clu, cross.mat = cross.mat,
                      par.mat = par.mat)
      
      # apply a constraint (remove 1st column) if par or anc Q.eff
      
      if((Q.eff == "par") | (Q.eff == "anc")) Q[, -1, drop = FALSE] else Q
      
      
    }
    
  }
  
  Q.list <- mapply(FUN = IncMat_QTL_MQE, x = Q.pos, Q.eff = Q.eff,
                   MoreArgs = list(mppData = mppData, mppData_bi = mppData_bi,
                                   par.clu = par.clu, cross.mat = cross.mat,
                                   par.mat = parent.mat), SIMPLIFY = FALSE)
  
  names(Q.list) <- paste0("Q", 1:length(Q.list))
  
  
  # 3. Compute the models
  #######################
  
  ind <- TRUE
  
  while(ind) {
    
    ### 3.1 elaboration of model formulas
    
    model.formulas <- formula_backward(Q.names = names(Q.list), VCOV = VCOV)
    
    ### 3.2 computation of the models
    
    pvals <- lapply(X = model.formulas, FUN = QTLModelBack, mppData = mppData,
                    Q.list = Q.list, cross.mat = cross.mat, VCOV = VCOV)
    
    
    pvals <- unlist(pvals)
    
    
    ### 3.4 test the p-values
    
    
    if(sum(pvals > alpha) > 0) {
      
      # remove the QTL position from the list
      
      Q.list <- Q.list[!(pvals==max(pvals))]
      
      # test if there is no more positions
      
      if(length(Q.list) == 0){ind <- FALSE}
      
    } else {
      
      # stop the procedure
      
      ind <- FALSE
      
    }
    
  }
  
  # 4. return the final list of QTLs
  ##################################
  
  index.rem <- as.numeric(substr(names(Q.list), 2, nchar(names(Q.list))))
  
  QTL.rem <- QTL[index.rem] ; Qeff.rem <- Q.eff[index.rem]
  
  QTL <- mppData$map[mppData$map[, 1] %in% QTL.rem, ]
  
  QTL <- data.frame(QTL, Qeff.rem, stringsAsFactors = FALSE)
  
  colnames(QTL)[5] <- "QTL.eff"
  
  if(dim(QTL)[1] == 0){
    
    QTL <- NULL
    
  }
  
  return(QTL)
  
}