##########
# QTL_R2 #
##########

#' QTL global and partial R squared
#' 
#' Computes the global and partial (adjusted) R squared of a  list of QTLs.
#' 
#' The function calculates two type of R squared. For all models
#' (\code{VCOV = "h.err", "h.err.as", "pedigree", "ped_cr.err"}), except
#' the CSRT model, the function computes the R squared using a linear model.
#' For the CSRT model \code{VCOV = "cr.err"}, it is possible to estimate a
#' likelihood ratio R squared as derived by Cox and Snell (1989).
#' 
#' In both situations, we estimate the extra variance explained by a full model
#' containing the QTL terms with respect to a reduced model containing only
#' the cross intercept terms. For the linear model R squared the function use
#' the ratio between the
#' residual sum of square of these two models: R2 = 1-(RSS(f))/(RSS(r)).
#' For the likelihood ratio R squared the full and reduced models are computed
#' using maximum likelihood and not REML following the recommentation of
#' Sun et al. (2010).
#' 
#' Partial R squared for each individual QTL position can also be calculated.
#' Two types of partial R squared are returned. The first one
#' uses the difference between the R squared obtained with all QTL
#' positions and the R squared obtain with all position minus the ith one
#' (difference R squared). The second method used only the ith QTL position
#' in the model (single R squared).
#' 
#' For both global and partial R squared, it is possible to obtained adjusted
#' measurements taking the number of degrees of freedom into consideration using
#' an adaptation of the formula given by Utz et al. (2000):
#' R.adj = R-(z/(N-z-n.cr))*(1-R) where z is the total
#' number of estimated components of the genetic effect. N is the total number
#' of phenotypic information, and n.cr is the number of intercept (cross) terms.
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
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker or inbetween marker positions names.
#' Default = NULL.
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
#' @param LR.R2 \code{Logical} value. If LR.R2 = TRUE, the likelihood ratio
#' R squared will be used if \code{VCOV = "cr.err"}. In all other situations
#' R squared from a linear model will be computed. Default = TRUE.
#' 
#' @param glb.only \code{Logical} value. If LR.R2 = TRUE, only the global and 
#' global adjusted R squared will be returned. Default = FALSE.
#'
#' @return Return:
#'
#' \code{List} containing the following objects:
#'
#' \item{glb.R2 }{ Global R squared of all QTL terms.}
#'
#' \item{glb.adj.R2 }{ Global adjusted R squared of all QTL terms.}
#'
#' \item{part.R2.diff }{ Vector of partial R squared doing
#' the difference between the full model and a model minus the ith QTL.}
#' 
#' \item{part.adj.R2.diff }{ Vector of partial adjusted R squared doing
#' the difference between the full model and a model minus the ith QTL.}
#' 
#' \item{part.R2.sg }{ Vector of partial R squared using only the ith QTL.}
#' 
#' \item{part.adj.R2.sg }{ Vector of partial adjusted R squared using only the
#' ith QTL.}
#' 
#' 
#' @author Vincent Garin
#'
#' @references
#' 
#' Cox, D. R., & Snell, E. J. (1989). Analysis of binary data (Vol. 32).
#' CRC Press.
#' 
#' Sun, G., Zhu, C., Kramer, M. H., Yang, S. S., Song, W., Piepho, H. P., & Yu, J.
#' (2010). Variation explained in mixed-model association mapping. Heredity,
#' 105(4), 333-340.
#' 
#' Utz, H. F., Melchinger, A. E., & Schon, C. C. (2000). Bias and sampling error
#' of the estimated proportion of genotypic variance explained by quantitative
#' trait loci determined from experimental data in maize using cross validation
#' and validation with independent samples. Genetics, 154(4), 1839-1849.
#'
#' @seealso \code{\link{mppData_form}}, \code{\link{parent_cluster}},
#' \code{\link{QTL_select}}, \code{\link{USNAM_parClu}}
#'
#' @examples
#'
#' data(USNAM_mppData)
#' 
#' SIM <- mpp_SIM(USNAM_mppData)
#' 
#' QTL <- QTL_select(Qprof = SIM, threshold = 3, window = 20)
#' 
#' QTL_R2(mppData = USNAM_mppData, QTL = QTL, Q.eff = "cr", VCOV = "h.err")
#' 
#' QTL_R2(mppData = USNAM_mppData, QTL = QTL, Q.eff = "cr", VCOV = "cr.err")
#' 
#' @export
#'

QTL_R2 <- function(mppData, QTL = NULL, Q.eff = "cr", par.clu = NULL,
                   VCOV = "h.err", LR.R2 = TRUE, glb.only = FALSE){
  
  # 1. check the data format
  ##########################
  
  check.model.comp(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
                   par.clu = par.clu, QTL = QTL, fct = "R2")
  
  if((VCOV == "cr.err") & (!LR.R2)) VCOV <-  "h.err"
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 inverse of the pedigree matrix
  
  # formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.2 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.3 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.4 modify the par.clu object order parents columns and replace
  # monomorphic
  
  if (Q.eff == "anc") {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.5 Formation of the list of QTL
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    
  }
  
  Q.list <- lapply(X = Q.pos, FUN = IncMat_QTL, mppData = mppData,
                   cross.mat = cross.mat, par.mat = parent.mat,
                   par.clu = par.clu, Q.eff = Q.eff)
  
  # re-arrange the QTL incidences matrices if parental or ancestral model
  
  if((Q.eff == "par") || (Q.eff == "anc")){
    Q.list <- lapply(X = Q.list, FUN = function(x) x[, -1, drop = FALSE])
  }
  
  # count the number of QTL elements for degree of freedom
  
  zi <- vapply(X = Q.list, function(x) dim(x)[2], FUN.VALUE = numeric(1))  
  
  n.QTL <- length(Q.list)
  
  # 3. Compute the R squared
  ##########################
  
  if(VCOV == "cr.err"){
    
    ### 3.1 Global adjusted and unadjusted linear R squared
    
    R2.all <- R2_LR(mppData = mppData, QTL = do.call(cbind, Q.list))
    R2 <- R2.all[1]
    R2.adj <- R2.all[2]
    
    if(glb.only) {
      
      return(list(glb.R2 = R2, glb.adj.R2 = R2.adj))
      
    } else {
      
      if(n.QTL > 1){
        
        ### 3.2 Compute the partial R squared
        
        # functions to compute the R squared or all QTL minus 1 or only 1 QTL position
        
        part.R2.diff <- function(x, QTL, mppData) {
          R2_LR(mppData = mppData, QTL = do.call(cbind, Q.list[-x]))[1]
        }
        
        part.R2.sg <- function(x, QTL, mppData) {
          R2_LR(mppData = mppData, QTL = do.call(cbind, Q.list[x]))[1]
        }
        
        R2_i.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff, QTL = Q.list,
                           mppData = mppData)
        
        R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
        
        R2_i.sg <- lapply(X = 1:n.QTL, FUN = part.R2.sg, QTL = Q.list,
                          mppData = mppData)
        
        R2_i.sg <- unlist(R2_i.sg)
        
        N <- sum(!is.na(mppData$trait[, 1]))
        
        # adjust the values
        
        
        R2_i.dif.adj <- R2_i.dif - ((zi/(N - zi - mppData$n.cr)) * (100 - R2_i.dif))
        
        R2_i.sg.adj <- R2_i.sg - ((zi/(N-zi-mppData$n.cr)) * (100 - R2_i.sg))
        
        names(R2_i.dif) <- names(R2_i.dif.adj) <- paste0("Q", 1:n.QTL)
        names(R2_i.sg) <- names(R2_i.sg.adj) <- paste0("Q", 1:n.QTL)
        
        return(list(glb.R2 = R2,
                    glb.adj.R2 = R2.adj,
                    part.R2.diff = R2_i.dif,
                    part.adj.R2.diff = R2_i.dif.adj,
                    part.R2.sg = R2_i.sg,
                    part.adj.R2.sg = R2_i.sg.adj))
        
      } else {
        
        names(R2) <- names(R2.adj) <- "Q1"
        
        return(list(glb.R2 = R2,
                    glb.adj.R2 = R2.adj,
                    part.R2.diff = R2,
                    part.adj.R2.diff = R2.adj,
                    part.R2.sg = R2,
                    part.adj.R2.sg = R2.adj))
      }
        
    }
    
    
    
  } else {
    
    # linear model and other VCOV (pedigree, pedigree_cr.err)
    
    ### 3.1 Global adjusted and unadjusted linear R squared
    
    R2.all <- R2_lin(mppData = mppData, QTL = do.call(cbind, Q.list))
    R2 <- R2.all[1]
    R2.adj <- R2.all[2]
    
    if(glb.only) {
      
      return(list(glb.R2 = R2, glb.adj.R2 = R2.adj))
      
    } else {
      
      if(n.QTL > 1){
        
        ### 3.2 Compute the partial R squared
        
        # functions to compute the R squared or all QTL minus 1 or only 1 QTL position
        
        part.R2.diff <- function(x, QTL, mppData) {
          R2_lin(mppData = mppData, QTL = do.call(cbind, Q.list[-x]))[1]
        }
        
        part.R2.sg <- function(x, QTL, mppData) {
          R2_lin(mppData = mppData, QTL = do.call(cbind, Q.list[x]))[1]
        }
        
        R2_i.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff, QTL = Q.list,
                           mppData = mppData)
        
        R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
        
        R2_i.sg <- lapply(X = 1:n.QTL, FUN = part.R2.sg, QTL = Q.list,
                          mppData = mppData)
        
        R2_i.sg <- unlist(R2_i.sg)
        
        N <- sum(!is.na(mppData$trait[, 1]))
        
        # adjust the values
        
        
        R2_i.dif.adj <- R2_i.dif - ((zi/(N - zi - mppData$n.cr)) * (100 - R2_i.dif))
        
        R2_i.sg.adj <- R2_i.sg - ((zi/(N-zi-mppData$n.cr)) * (100 - R2_i.sg))
        
        names(R2_i.dif) <- names(R2_i.dif.adj) <- paste0("Q", 1:n.QTL)
        names(R2_i.sg) <- names(R2_i.sg.adj) <- paste0("Q", 1:n.QTL)
        
        return(list(glb.R2 = R2,
                    glb.adj.R2 = R2.adj,
                    part.R2.diff = R2_i.dif,
                    part.adj.R2.diff = R2_i.dif.adj,
                    part.R2.sg = R2_i.sg,
                    part.adj.R2.sg = R2_i.sg.adj))
        
      } else {
        
        names(R2) <- names(R2.adj) <- "Q1"
        
        return(list(glb.R2 = R2,
                    glb.adj.R2 = R2.adj,
                    part.R2.diff = R2,
                    part.adj.R2.diff = R2.adj,
                    part.R2.sg = R2,
                    part.adj.R2.sg = R2.adj))
      }
      
    }
    
  } 
    
}