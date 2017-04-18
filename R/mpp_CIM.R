###########
# mpp_CIM #
###########

#' MPP composite interval maping
#' 
#' Computes QTL models using cofactors with different possible assumptions
#' concerning the number of alleles at the QTL position and the variance
#' covariance structure (VCOV) of the model. For more details about the
#' different models, see documentation of the function \code{\link{mpp_SIM}}.
#' The function returns a -log10(p-value) QTL profile.
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
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. For more details see \code{\link{mpp_SIM}}. Default = "cr".
#'
#' @param par.clu Required argument for the ancesral model \code{(Q.eff = "anc")}.
#' \code{interger matrix} representing the results of a parents genotypes
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
#' @param cofactors Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker or inbetween marker positions names.
#' Default = NULL.
#' 
#' @param window \code{Numeric} distance on the left an right of a cofactor
#' position where it is not included in the model. Default = 20.
#' 
#' @param est.gen.eff \code{Logical} value. if \code{est.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot_genEffects}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' \strong{This functionality is ony available for the cross-specific,
#' parental and ancestral models.}
#' Default value = FALSE.
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
#' \item{CIM }{\code{Data.frame} of class \code{QTLprof}. with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) Interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-val). And if
#' \code{est.gen.eff = TRUE}, p-values of the cross or parental QTL effects.}
#' 
#' @author Vincent Garin
#' 
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mpp_SIM}},
#' \code{\link{parent_cluster}}, \code{\link{QTL_select}},
#' \code{\link{USNAM_parClu}}
#' 
#' 
#' @examples
#' 
#' # data
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' 
#' # Cross-specific effect model
#' #############################
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData, Q.eff = "cr", VCOV = "h.err")
#' 
#' cofactors <- QTL_select(Qprof = SIM, threshold = 3, window = 20)
#' 
#' CIM <- mpp_CIM(mppData = USNAM_mppData, Q.eff = "cr", VCOV = "h.err",
#' cofactors = cofactors, window = 20, est.gen.eff = TRUE)
#' 
#' plot_QTLprof(Qprof = CIM)
#' plot_genEffects(mppData = USNAM_mppData, Qprof = CIM, Q.eff = "cr")
#' 
#' \dontrun{
#' 
#' # Using multiple core
#' 
#' library(parallel)
#' n.cores <- detectCores()
#' cluster <- makeCluster((n.cores-1))
#' 
#' CIM <- mpp_CIM(mppData = USNAM_mppData, Q.eff = "cr", VCOV = "h.err",
#' cofactors = cofactors, window = 20, est.gen.eff = TRUE, parallel = TRUE,
#' cluster = cluster)
#' 
#' stopCluster(cl = cluster)
#' 
#' }
#' 
#' # Bi-allelic model
#' ##################
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData_bi, Q.eff = "biall", VCOV = "h.err")
#' 
#' cofactors <- USNAM_mppData_bi$map[c(15, 63), 1]
#' 
#' CIM <- mpp_CIM(mppData = USNAM_mppData_bi, Q.eff = "biall", VCOV = "h.err",
#' cofactors = cofactors, window = 20)
#' 
#' plot_QTLprof(Qprof = CIM, type = "h")
#'                                
#' @export
#' 


mpp_CIM <- function(mppData, Q.eff = "cr", par.clu = NULL, VCOV = "h.err",
                    cofactors = NULL,  window = 20, est.gen.eff = FALSE,
                    parallel = FALSE, cluster = NULL)
{
  
  # 1. Check data format and arguments
  ####################################
  
  check.model.comp(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
                   par.clu = par.clu, est.gen.eff = est.gen.eff,
                   parallel = parallel, cluster = cluster,
                   cofactors = cofactors, fct = "CIM")
  
  # 2. Form required elements for the analysis
  ############################################
  
  ### 2.1 inverse of the pedigree matrix
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.2 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.3 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.4 modify the par.clu object order parents columns and replace monomorphic
  
  if (Q.eff == "anc") {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.5 Formation of the list of cofactors
  
  if(is.character(cofactors)){
    
    cof.pos <- which(mppData$map[, 1] %in% cofactors)
    
  } else {
    
    cof.pos <- which(mppData$map[, 1] %in% cofactors[, 1])
    
  }
  
  cof.list <- lapply(X = cof.pos, FUN = IncMat_QTL, mppData = mppData,
                     cross.mat = cross.mat, par.mat = parent.mat,
                     par.clu = par.clu, Q.eff = Q.eff, order.MAF = TRUE)
  
  ### 2.6 Formation of the genome-wide and cofactors partition
  
  vect.pos <- 1:dim(mppData$map)[1]
  
  # 2.6.1 cofactor partition tells if the cofactor should be included or
  # not in the model at each position.
  
  if (is.character(cofactors)){
    
    cofactors2 <- mppData$map[mppData$map[, 1] %in% cofactors, c(2, 4)]
    
  } else { cofactors2 <- cofactors[, c(2, 4)] }
  
  test.cof <- function(x, map, window) {
    
    t1 <- map$chr == as.numeric(x[1])
    t2 <- abs(map$pos.cM - as.numeric(x[2])) < window
    !(t1 & t2)
    
  }
  
  cof.part <- apply(X = cofactors2, MARGIN = 1, FUN = test.cof,
                    map = mppData$map, window = window)
  
  
  # 3. computation of the CIM profile (genome scan)
  #################################################
  
  if (parallel) {
    
    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelCIM,
                          mppData = mppData, cross.mat = cross.mat,
                          par.mat = parent.mat, Q.eff = Q.eff,
                          par.clu = par.clu, VCOV = VCOV, cof.list = cof.list,
                          cof.part = cof.part, est.gen.eff = est.gen.eff)
    
  } else {
    
    log.pval <- lapply(X = vect.pos, FUN = QTLModelCIM,
                       mppData = mppData, cross.mat = cross.mat,
                       par.mat = parent.mat, Q.eff = Q.eff,
                       par.clu = par.clu, VCOV = VCOV, cof.list = cof.list,
                       cof.part = cof.part, est.gen.eff = est.gen.eff)
    
  }
  
  
  log.pval <- t(data.frame(log.pval))
  if(est.gen.eff & (VCOV == "h.err")){log.pval[is.na(log.pval)] <- 1}
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0
  
  # 4. form the results
  #####################
  
  CIM <- data.frame(mppData$map, log.pval)
  
  
  if(est.gen.eff){
    
    if(Q.eff == "cr"){ Qeff_names <- unique(mppData$cross.ind)
    
    } else { Qeff_names <- mppData$parents }
    
    colnames(CIM)[5:dim(CIM)[2]] <- c("log10pval", Qeff_names)
    
  } else {colnames(CIM)[5] <- "log10pval"}
  
  class(CIM) <- c("QTLprof", "data.frame")
  
  ### 4.1: Verify the positions for which model could not be computed
  
  if(sum(CIM$log10pval == 0) > 0) {
    
    if (sum(CIM$log10pval) == 0){
      
      message("The computation of the QTL models failled for all positions.
              This could be due to problem in asreml() function.")
      
    } else {
      
      list.pos <- mppData$map[(CIM$log10pval == 0), 1]
      
      end.mess <- ". This could be due to singularities or function issue."
      
      text <- paste("The computation of the QTL model failed for the following positions: ",
                    paste(list.pos, collapse = ", "), end.mess)
      
      message(text)
      
    }
    
  }
  
  return(CIM)
  
}