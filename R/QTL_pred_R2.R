###############
# QTL_pred_R2 #
###############

#' Predicted QTL global and partial R squared
#' 
#' Computes predicted R squared in a validation set using QTLs detected in a
#' training set.
#'
#' Compute QTLs predicted R squared in a validation set (\code{mppData.vs}).
#' These QTLs have been previously detected in a training set
#' (\code{mppData.ts}). The global R squared is obtained using the Pearson
#' squared correlation between the observed trait values in the validation set
#' (y.vs) and predicted values using estimated QTL effects in the training set
#' (y.pred.vs = X.vs * B.ts). The correlation can be calculated within cross and
#' then averaged (\code{within.cross = TRUE}) or at the whole population level.
#' 
#' Partial R squared statistics are also calculated for each individual position.
#' Two types or partial R squared are computed. The first one making the
#' difference between the global R squared and the R squared computed without
#' the ith position (difference R squared). The second method only uses
#' the ith QTL for trait values predition (single R squared) (for details
#' see documentation of the function \code{\link{QTL_R2}}).
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
#' @param mppData.ts An object of class \code{mppData} for the training set. See
#' \code{\link{mppData_form}} for details.
#' 
#' @param mppData.vs An object of class \code{mppData} for the validation set.
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
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker or in between marker positions names.
#' Default = NULL.
#' 
#' @param within.cross \code{Logical} value indicating if the predicted
#' R squared must be computed at the whole population level or within cross.
#' In the later case, predicted R squared are computed within cross and the
#' average is returned. Default = TRUE.
#' 
#' 
#' @return Return:
#' 
#' \code{List} containing the following objects:
#'
#' \item{glb.R2 }{ Global predicted R squared of all QTL terms.}
#'
#' \item{part.R2.diff }{ Vector of predicted partial R squared doing
#' the difference between the full model and a model minus the ith QTL.}
#' 
#' \item{part.R2.sg }{ Vector of predicted partial R squared using only the
#' ith QTL.}
#' 
#' 
#' @author Vincent Garin
#' 
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{parent_cluster}},
#' \code{\link{QTL_R2}}, \code{\link{QTL_select}}, \code{\link{USNAM_parClu}}
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' 
#' folds <- CV_partition(cross.ind = USNAM_mppData$cross.ind, k = 5)
#' 
#' mppData.ts <- mppData_subset(mppData = USNAM_mppData,
#'                              gen.list = folds[[1]]$train.set)
#' 
#' mppData.vs <- mppData_subset(mppData = USNAM_mppData,
#'                              gen.list = folds[[1]]$val.set)
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData)
#' QTL <- QTL_select(SIM)
#'
#' QTL_pred_R2(mppData.ts = mppData.ts, mppData.vs = mppData.vs, QTL = QTL)
#' 
#' @export
#' 


QTL_pred_R2 <- function(mppData.ts, mppData.vs, Q.eff = "cr", par.clu = NULL,
                        VCOV = "h.err", QTL = NULL, within.cross = TRUE) {
  
  # 1. test data format
  #####################
  
  check.model.comp(Q.eff = Q.eff, VCOV = VCOV, par.clu = par.clu, QTL = QTL,
                   mppData.ts = mppData.ts, mppData.vs = mppData.vs,
                   fct = "R2_pred")
  
  
  if(is.character(QTL)){ n.QTL <- length(QTL) } else { n.QTL <- dim(QTL)[1] }
  
  # 2. obtain the genetic effects (Betas)
  #######################################
  
  if((Q.eff == "biall") || (Q.eff == "cr")){
    if(!is.null(mppData.ts$geno.par)) {mppData.ts$geno.par <- NULL}
  }
  
  effects <- QTL_genEffects(mppData = mppData.ts, QTL = QTL, Q.eff = Q.eff,
                            par.clu = par.clu, VCOV = VCOV)
  
  # need to re-order the row of the effects according to the parent list
  
  if(Q.eff == "par"){
    
    eff.names <- substr(rownames(effects[[1]]), 3 , nchar(rownames(effects[[1]])))
    index <- match(eff.names, mppData.vs$parents)
    
    effects <- lapply(X = effects, function(x, ind) x[ind, ], ind = index)
    
  } else if (Q.eff == "anc"){
    
    effects <- lapply(X = effects, function(x, ind) x[ind, ],
                      ind = mppData.vs$parents)
    
  }
  
  
  B.ts <- lapply(X = seq_along(effects), FUN = function(x, Qeff) Qeff[[x]][, 1],
                 Qeff = effects)
  
  # 3. obtain the QTL incidence matrices of the positions (X.vs)
  ##############################################################
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData.vs$map[, 1] %in% QTL)
    
  } else {
    
    Q.pos <- which(mppData.vs$map[, 1] %in% QTL[, 1])
    
  }
  
  # switch Qeff for parental QTL incidence matrix if ancestral model
  
  if(Q.eff == "anc") Q.eff.part <- "par" else Q.eff.part <- Q.eff
  
  Q.list <- lapply(X = Q.pos, FUN = IncMat_QTL, mppData = mppData.vs,
                   cross.mat = IncMat_cross(mppData.vs$cross.ind),
                   par.mat = IncMat_parent(mppData.vs), par.clu = par.clu,
                   Q.eff = Q.eff.part)
  
  # 4. Predicted R squared computation cor(y.vs, X.vs * B.ts)^2
  ##############################################################
  
  # global R squared
  
  R2 <- R2_pred(mppData.vs = mppData.vs, B.ts = B.ts, Q.list = Q.list,
                within.cross = within.cross)
  
  # partial R2
  
  if (n.QTL > 1) {
    
    part.R2.diff <- function(x, mppData.vs, B.ts, Q.list, within.cross) {
      R2_pred(mppData.vs = mppData.vs, B.ts = B.ts[-x], Q.list =  Q.list[-x],
              within.cross = within.cross)
    }
    
    part.R2.sg <- function(x, mppData.vs, B.ts, Q.list, within.cross) {
      R2_pred(mppData.vs = mppData.vs, B.ts = B.ts[x], Q.list =  Q.list[x],
              within.cross = within.cross)
    }
    
    R2_i.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff, mppData.vs = mppData.vs,
                       B.ts = B.ts, Q.list = Q.list,
                       within.cross = within.cross)
    
    R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
    
    R2_i.sg <- lapply(X = 1:n.QTL, FUN = part.R2.sg, mppData.vs = mppData.vs,
                      B.ts = B.ts, Q.list = Q.list,
                      within.cross = within.cross)
    
    R2_i.sg <- unlist(R2_i.sg)
    
    
    names(R2_i.dif) <- paste0("Q", 1:n.QTL)
    names(R2_i.sg) <- paste0("Q", 1:n.QTL)
    
    
    return(list(glb.R2 = R2, part.R2.diff = R2_i.dif, part.R2.sg = R2_i.sg))
    
  } else {
    
    names(R2) <- "Q1"
    
    return(list(glb.R2 = R2, part.R2.diff = R2, part.R2.sg = R2))
    
  }
  
}