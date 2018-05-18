################
# mpp_back_elim #
################

#' Backward elimination on QTL candidates
#' 
#' Performs a backward elimination using a list of given QTLs positions. The
#' positions with a p-value above the significance level \code{alpha}, are
#' successively removed.
#' 
#' The function starts with all QTL positions in the model and test the inclusion
#' of each position as the last in the model. If all position p-values are below
#' \code{alpha} the procedure stop. If not the position with the highest p-value
#' is remove and the procedure continue until there is no more unsignificant
#' position.
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
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#'
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker or in between marker positions names.
#' Default = NULL.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
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
#' \item{QTL }{\code{Data.frame} of class \code{QTLlist} with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-values).}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mpp_SIM}},
#' \code{\link{parent_cluster}}, \code{\link{USNAM_parClu}}
#' 
#' @examples
#' 
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData)
#' 
#' QTL <- QTL_select(SIM)
#' 
#' QTL.sel <- mpp_back_elim(mppData = mppData, QTL = QTL)
#' 
#' @export
#' 


mpp_back_elim <- function (mppData, trait = 1, QTL = NULL, Q.eff = "cr",
                         VCOV = "h.err", alpha = 0.05) {
  
  # 1. Check data format
  ######################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                   QTL = QTL, fct = "back")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 trait values
  
  if(is.numeric(trait)){
    
    t_val <- mppData$pheno[, trait]
    
  } else {
    
    trait.names <- colnames(mppData$pheno)
    t_val <- mppData$pheno[, which(trait %in% trait.names)]
    
  }
  
  ### 2.2 inverse of the pedigree matrix
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 Formation of the list of QTL
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    
    QTL <- mppData$map[mppData$map[, 1] %in% QTL, ]
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    
  }
  
  Q.list <- lapply(X = Q.pos, FUN = inc_mat_QTL, mppData = mppData,
                    Q.eff = Q.eff, order.MAF = TRUE)
  
  names(Q.list) <- paste0("Q", 1:length(Q.list))
  
  
  # 3. Compute the models
  #######################
  
  ind <- TRUE
  
  while(ind) {
    
    ### 3.1 elaboration of model formulas
    
    model.formulas <- formula_backward(Q.names = names(Q.list), VCOV = VCOV)
    
    ### 3.2 computation of the models
    
    pvals <- lapply(X = model.formulas, FUN = QTLModelBack, mppData = mppData,
                    trait = t_val, Q.list = Q.list, cross.mat = cross.mat,
                    VCOV = VCOV)
    
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
  
  QTL <- QTL[as.numeric(substr(names(Q.list), 2, nchar(names(Q.list)))), ,
             drop = FALSE]
  
  if(dim(QTL)[1] == 0){
    
    QTL <- NULL
    
    return(QTL)
    
  } else{
    
    class(QTL) <- c("QTLlist", "data.frame")
    
    return(QTL)
    
  }
  
}