##################
# QTL_genEffects #
##################

#' QTL genetic effects
#' 
#' Computes a multi-QTL model with a list QTL candidates(\code{QTL}) and return
#' the decomposed QTL effects per cross or per parents.
#' 
#' This function computes for each QTL position the genetic effects of the
#' cross, parental, ancestral or SNP allele components. For the cross-specific
#' model (\code{Q.eff = "cr"}), the genetics effects represent the substitution
#' effect of an single allele from the parent 2 (or B) with respect to an allele
#' coming from the parent 1 or A. All effects are given in absolute value with
#' the parent that cary the positive allele.
#' 
#' For the parental and the ancestral model (\code{Q.eff = "par" or "anc"}), it
#' is possible to estimate maximum n-1 parental or ancestral alleles per
#' interconnected part of the design. For these two models, one
#' parental (ancestral) allele is set as reference per interconnected part of the
#' design. Effects of the other alleles are estimated as deviation with respect
#' to the reference. Connected parts of the design can be determined using Weeks
#' and Williams (1964) method (\code{\link{design_connectivity}}). By default,
#' the reference allele is the most frequent one. The user can also specify a
#' parental allele that will be used as reference using the argument
#' \code{ref.par}. This option is only available if the MPP design is composed
#' of a unique connected part.
#' 
#' For the bi-allelic model (\code{Q.eff = "biall"}), the genetic effects
#' represent the effects of a single allele copy of the most frequent allele.
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
#' @param ref.par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. For the
#' ancestral model, the ancestral class containing the reference parent will be
#' set as reference. \strong{This option can only be used if the MPP design is
#' composed of a unique connected part}. Default = NULL.
#' 
#' 
#' @return Return:
#'
#' \item{Qeff}{\code{List} of \code{data.frame} (one per QTL) containing the
#' following information:
#' 
#' \enumerate{
#' 
#' \item{QTL genetic effects per cross or parent.}
#' \item{Standard error of the QTL effects.}
#' \item{Test statistics of the effects (t-test or Wald statistic).}
#' \item{P-value of the test statistics.}
#' \item{Significance of the QTL effects.}
#' \item{For cross-specific model, parent with the positive additive effects.}
#' \item{For parental and ancestral model, indicator of connected part of the
#' design and reference.}
#' \item{Allele scores of the parents if \code{geno.par} is non NULL
#' in the \code{mppData} object.}
#' 
#' }
#' 
#' }
#' 
#' \item{tab.Qeff}{\code{data.frame} with one column per QTL giving the
#' QTL genetic effect per cross or per parent with its significance. The
#' first two rows indicate the chromosome and the position in cM of each
#' QTL.}
#' 
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{parent_cluster}},
#' \code{\link{QTL_select}}, \code{\link{USNAM_parClu}}
#' 
#' @references 
#' 
#' Weeks, D. L., & Williams, D. R. (1964). A note on the determination of
#' connectedness in an N-way cross classification. Technometrics, 6(3), 319-324.
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' data(USNAM_parClu)
#' par.clu <- USNAM_parClu
#' 
#' # QTL candidates
#' 
#' SIM <- mpp_SIM(USNAM_mppData)
#' QTL <- QTL_select(SIM)
#' 
#' # cross-specific model
#' 
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData, QTL = QTL, Q.eff = "cr")
#' QTL.effects
#' 
#' # parental model
#' 
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData, QTL = QTL,
#'                                Q.eff = "par")
#' QTL.effects
#' 
#' # ancestral model
#' 
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData, QTL = QTL, Q.eff = "anc",
#'                               par.clu = par.clu)
#' QTL.effects
#' 
#' # bi-allelic model
#' 
#' SIM <- mpp_SIM(USNAM_mppData_bi, Q.eff = "biall")
#' QTL <- QTL_select(SIM)
#' 
#' 
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData_bi, QTL = QTL,
#' Q.eff = "biall")
#' 
#' QTL.effects
#' 
#' @export 
#'


QTL_genEffects <- function(mppData, QTL = NULL, Q.eff = "cr", par.clu = NULL,
                           VCOV = "h.err", ref.par = NULL) {
  
  # 1. Check data format
  ######################
  
  check.model.comp(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
                   par.clu = par.clu, QTL = QTL, ref.par = ref.par,
                   fct = "QTLeffects")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 Phenotypic values
  
  trait <- mppData$trait[, 1]
  
  ### 2.2 inverse of the pedigree matrix
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.5 modify the par.clu object order parents columns and replace monomorphic
  
  if (Q.eff == "anc") {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.6 Formation of the list of QTL
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    
  }
  
  Q.list <- lapply(X = Q.pos, FUN = IncMat_QTL, mppData = mppData,
                   cross.mat = cross.mat, par.mat = parent.mat,
                   par.clu = par.clu, Q.eff = Q.eff)
  
  names(Q.list) <- paste0("Q", 1:length(Q.list))
  
  
  ### 2.7 For the parental and ancestral model organise the matrix to get the
  # desired constraint.
  
  Q.eff_temp <- rep(Q.eff, length(Q.list))
  
  order.Qmat <- mapply(FUN = IncMat_QTL_MAF, QTL = Q.list,
                       Q.eff_i = Q.eff_temp, Q.pos_i = Q.pos,
                       MoreArgs = list(mppData = mppData, par.clu = par.clu,
                                       ref.par = ref.par),
                       SIMPLIFY = FALSE)
  
  
  Q.list <- lapply(X = order.Qmat, FUN = function(x) x$QTL)
  allele_order <- lapply(X = order.Qmat, FUN = function(x) x$allele_order)
  con.ind <- lapply(X = order.Qmat, FUN = function(x) x$con.ind)
  
  
  # 3. model computation
  ######################
  
  model <- QTLModelQeff(mppData = mppData, trait = trait, cross.mat = cross.mat,
                        Q.list = Q.list, VCOV = VCOV)
  
  
  # 4. data processing
  ####################
  
  
  results <- Qeff_res_processing(model = model, mppData = mppData,
                                 cross.mat =  cross.mat, Q.list = Q.list,
                                 QTL = QTL, Q.eff = Q.eff, par.clu = par.clu,
                                 VCOV = VCOV, allele_order = allele_order,
                                 con.ind = con.ind)
  
  
  names(results) <- paste0("Q", 1:length(results))
  
  # table with all QTL effect per cross or per parents.
  
  if (Q.eff == "cr"){
    
    Qeff_sign <- lapply(results, `[`, c(1, 5))
    Qeff_sign <- lapply(Qeff_sign, function(x) paste(round(x[, 1], 3), x[, 2]))
    table.QTL <- data.frame(Qeff_sign)
    pos.info <- t(data.frame(QTL[, 2], QTL[, 4]))
    colnames(pos.info) <- paste0("Q", 1:length(results))
    table.QTL <- rbind.data.frame(pos.info, table.QTL)
    rownames(table.QTL) <- c("chr", "pos.cM", unique(mppData$cross.ind))
    
  } else if (Q.eff == "biall"){
    
    Qeff_sign <- lapply(results, `[`, c(1, 5))
    Qeff_sign <- lapply(Qeff_sign, function(x) paste(round(x[, 1], 3), x[, 2]))
    table.QTL <- data.frame(Qeff_sign)
    pos.info <- t(data.frame(QTL[, 2], QTL[, 4]))
    colnames(pos.info) <- paste0("Q", 1:length(results))
    table.QTL <- rbind.data.frame(pos.info, table.QTL)
    
    if (is.null(mppData$geno.par)){
      
      rownames(table.QTL) <- c("chr", "pos.cM", "Q.eff")
      
    } else {
      
      rownames(table.QTL) <- c("chr", "pos.cM", mppData$parents)
      
    }
    
  } else { # parental or ancestral
    
    Qeff_sign <- lapply(results, `[`, c(1, 5))
    
    # order according to parents
    
    Qeff_sign <- lapply(X = Qeff_sign, function(x, ind) x[ind, ],
                        ind = mppData$parents)
    
    Qeff_sign <- lapply(Qeff_sign, function(x) paste(round(x[, 1], 3), x[, 2]))
    table.QTL <- data.frame(Qeff_sign)
    pos.info <- t(data.frame(QTL[, 2], QTL[, 4]))
    colnames(pos.info) <- paste0("Q", 1:length(results))
    table.QTL <- rbind.data.frame(pos.info, table.QTL)
    rownames(table.QTL) <- c("chr", "pos.cM", mppData$parents)
    
  }
  
  return(list(Qeff = results, tab.Qeff = table.QTL))
  
}