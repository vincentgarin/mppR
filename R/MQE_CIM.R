###########
# MQE_CIM #
###########

#' Composite interval maping for MQE model
#' 
#' Computes a QTL profile with cofactors position having different type of QTL
#' effects.
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
#' @param Q.eff \code{Character} expression indicating type of QTL effect at the
#' tested position. Possibility to choose among: "cr", "par", "anc" or "biall".
#' For details look at \code{\link{mpp_SIM}}. Default = "cr".
#' 
#' @param par.clu Required argument if the user wants to allow QTLs with an
#' ancestral effect. \code{interger matrix} representing the results of a parents genotypes
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
#' @param cofactors Vector of \code{character} marker or in between marker
#' positions names. Default = NULL.
#' 
#' @param cof.Qeff \code{Character} vector indicating for each cofactor position
#' the type of QTL effect from: "cr", "par", "anc" and "biall".
#' 
#' @param chg.Qeff \code{Logical} value. If \code{chg.Qeff = TRUE}.
#' The type of QTL effect of the tested position will change for the one
#' specified in cof.Qeff when the function enter the region of a cofactor
#' delimited by \code{window}. Default = FALSE.
#' 
#' @param window \code{Numeric} value in centi-Morgan representing the distance
#' on the left an right of a cofactor position where it is not included in the
#' model. Default value = 20.
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
#'   
#' @return Return:
#' 
#' \item{CIM }{\code{Data.frame} of class \code{QTLprof} with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) Interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-values)}
#' 
#' @author Vincent Garin
#' 
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mpp_SIM}},
#' \code{\link{parent_cluster}}, \code{\link{QTL_select}}
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
#' # Equalize the list of markers of the two mppData objects and par.clu
#' 
#' com.mk.list <- intersect(mppData$map$mk.names, mppData_bi$map$mk.names)
#' 
#' mppData <- mppData_subset(mppData = mppData, mk.list = com.mk.list)
#' mppData_bi <- mppData_subset(mppData = mppData_bi, mk.list = com.mk.list)
#' par.clu <- par.clu[rownames(par.clu) %in% com.mk.list, ]
#' 
#' SIM <- mpp_SIM(mppData = mppData)
#' cofactors <- QTL_select(SIM)[, 1]
#' 
#' CIM <- MQE_CIM(mppData = mppData, mppData_bi = mppData_bi, Q.eff = "cr",
#'                par.clu = par.clu, cofactors = cofactors,
#'                cof.Qeff = c("anc", "par", "biall"))
#'
#'plot_QTLprof(CIM)
#'                                
#' @export
#' 



MQE_CIM <- function(mppData = NULL, mppData_bi = NULL, Q.eff = "cr",
                    par.clu = NULL, VCOV = "h.err", cofactors = NULL, cof.Qeff,
                    chg.Qeff = FALSE, window = 20, parallel = FALSE,
                    cluster = NULL){
  
  # 1. Check data format and arguments
  ####################################
  
  check.MQE(mppData = mppData, mppData_bi = mppData_bi, Q.eff = Q.eff,
            VCOV = VCOV, par.clu = par.clu, cofactors = cofactors,
            cof.Qeff = cof.Qeff, parallel = parallel, cluster = cluster,
            fct = "CIM")
  
  # 2. Form required elements for the analysis
  ############################################
  
  ### 2.1 inverse of the pedigree matrix
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.2 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.3 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.4 modify the par.clu object order parents columns and replace monomorphic
  
  if (("anc" %in% cof.Qeff)|| (Q.eff == "anc")) {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.5 Formation of the list of cofactors
  
  # order the list of cofactors and the corresponding Q.eff
  
  cof.pos <- vapply(X = cofactors,
                    FUN = function(x, mppData) which(mppData$map[,1] == x),
                    FUN.VALUE = numeric(1), mppData = mppData)
  
  cof.ord <- data.frame(cofactors, cof.Qeff, cof.pos,
                        stringsAsFactors = FALSE)
  
  cof.ord <- cof.ord[order(cof.pos),]
  
  cofactors <- cof.ord[, 1]; cof.Qeff <- cof.ord[, 2]
  
  cof.pos <- which(mppData$map[, 1] %in% cofactors)
  
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
  
  cof.list <- mapply(FUN = IncMat_QTL_MQE, x = cof.pos, Q.eff = cof.Qeff,
                     MoreArgs = list(mppData = mppData, mppData_bi = mppData_bi,
                                     par.clu = par.clu, cross.mat = cross.mat,
                                     par.mat = parent.mat), SIMPLIFY = FALSE)
  
  
  ### 2.6 Formation of the genome-wide and cofactors partition
  
  vect.pos <- 1:dim(mppData$map)[1]
  
  # 2.6.1 cofactor partition tells if the cofactor should be included or
  # not in the model at each position.
  
  cofactors2 <- mppData$map[cof.pos, c(2,4)]
  
  test.cof <- function(x, map, window) {
    
    t1 <- map$chr == as.numeric(x[1])
    t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
    !(t1 & t2)
    
  }
  
  cof.part <- apply(X = cofactors2, MARGIN = 1, FUN = test.cof,
                    map = mppData$map, window = window)
  
  # make a QTL effect partition to change type of QTL effect of the tested pos.
  
  if(chg.Qeff) {
    
    cof.part2 <- (!cof.part)*1
    
    Qeff.partition <- function(x, Q.eff, cof.Qeff){
      if(sum(x) == 0){Q.eff } else { cof.Qeff[max(which(x == 1))] }
    }
    
    Qeff.part <- apply(X = cof.part2, MARGIN = 1, FUN = Qeff.partition,
                       Q.eff = Q.eff, cof.Qeff = cof.Qeff)
    
  } else {Qeff.part <- rep(Q.eff, length(vect.pos))}
  
  
  # 3. computation of the CIM profile (genome scan)
  #################################################
  
  if (parallel) {
    
    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelCIM_MQE,
                          mppData = mppData, mppData_bi = mppData_bi, 
                          cross.mat = cross.mat, par.mat = parent.mat,
                          Qeff.part = Qeff.part, par.clu = par.clu, VCOV = VCOV,
                          cof.list = cof.list, cof.part = cof.part)
    
    
  } else {
    
    log.pval <- lapply(X = vect.pos, FUN = QTLModelCIM_MQE,
                       mppData = mppData, mppData_bi = mppData_bi, 
                       cross.mat = cross.mat, par.mat = parent.mat,
                       Qeff.part = Qeff.part, par.clu = par.clu, VCOV = VCOV,
                       cof.list = cof.list, cof.part = cof.part)
  }
  
  
  log.pval <- t(data.frame(log.pval))
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0
  
  # 4. form the results
  #####################
  
  CIM <- data.frame(mppData$map, log.pval)
  colnames(CIM)[5] <- "log10pval"
  class(CIM) <- c("QTLprof", "data.frame")
  
  ### 4.1: Verify the positions for which model could not be computed
  
  if(sum(CIM$log10pval == 0) > 0) {
    
    if (sum(CIM$log10pval) == 0){
      
      warning("The computation of the QTL models failled for all positions.
              This could be due to problem in asreml() function.")
      
    } else {
      
      list.pos <- mppData$map[(CIM$log10pval == 0), 1]
      
      end.mess <- ". This could be due to singularities or function issue."
      
      text <- paste("The computation of the QTL model failed for the following positions: ",
                    paste(list.pos, collapse = ", "), end.mess)
      
      warning(text)
      
    }
    
  }
  
  return(CIM)
  
}