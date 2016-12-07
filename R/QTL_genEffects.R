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
#' For the parental or ancestral models (\code{Q.eff = "par"
#' or "anc"}), the genetic effects represent the effect of a single parental
#' (ancestral) allele with respect to a reference allelele. The reference
#' parental (ancestral) allele can be specified in the argument
#' \code{par.ref}. For the ancestral model, the ancestral class inferred
#' for the specified parent will be used as reference. For the parental and
#' ancestral model computed with the homogeneous residual term assumption
#' (\code{VCOV = "h.err"}), it is also possible to estimate one effect for
#' each parental (ancestral) allele with the constraint that these effects
#' sum to zero. This can be done using (\code{const = "sum.0"}).
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
#' @param const \code{Character} expression indicating the type of contstraint
#' used for the estimation of the genetic effect of the HRT (linear) model
#' (\code{VCOV = "h.err"}) for the parental or ancestral models
#' (\code{Q.eff = "par" or "anc"}). If const = "set.0", the value of the
#' parent or the ancestral class containing the parent specified in argument
#' (\code{par.ref}) is set to 0. If const = "sum.0", the sum of
#' the parental (ancestral) effects is constrained to sum to 0.
#' Default = "set.0".
#' 
#' @param par.ref \code{Character} expression indicating the parent or the
#' ancestral class containing the specified parent that will be set as reference
#' for the computation of the genetic effect. By default
#' the function will use the first parent of the \code{mppData$parents} list.
#' 
#' 
#' @return Return:
#'
#' \item{results}{\code{List} of \code{data.frame} (one per QTL) containing the
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
#' \item{Allele scores of the parents if \code{geno.par} is non NULL
#' in the \code{mppData} object.}
#' 
#' }
#' 
#' }
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{parent_cluster}},
#' \code{\link{QTL_select}}, \code{\link{USNAM_parClu}}
#' 
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
#' # constraint(set one parents as reference - CML322)
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData, QTL = QTL, Q.eff = "par",
#'                               const = "set.0", par.ref = "CML322")
#' QTL.effects
#' 
#' # constraint(sum to zero)
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData, QTL = QTL, Q.eff = "par",
#'                               const = "sum.0")
#' QTL.effects
#' 
#' # ancestral model
#' 
#' # constraint(set one ancestral class as reference - M37W)
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData, QTL = QTL, Q.eff = "anc",
#'                               par.clu = par.clu, const = "set.0", par.ref = "M37W")
#' QTL.effects
#' 
#' # constraint(sum to zero)
#' QTL.effects <- QTL_genEffects(mppData = USNAM_mppData, QTL = QTL, Q.eff = "anc",
#'                               par.clu = par.clu, const = "sum.0")
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
                           VCOV = "h.err",  const = "set.0", par.ref = NULL) {
  
  # 1. Check data format
  ######################

  if(is.null(par.ref)){ par.ref <- mppData$parents[1]}

  check.model.comp(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
                   par.clu = par.clu, par.ref = par.ref, const = const,
                   QTL = QTL, fct = "QTLeffects")


  # 2. elements for the model
  ###########################

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


  ### 2.5 Formation of the list of QTL

  if(is.character(QTL)){

    Q.pos <- which(mppData$map[, 1] %in% QTL)

  } else {

    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])

  }

  Q.list <- lapply(X = Q.pos, FUN = IncMat_QTL, mppData = mppData,
                   cross.mat = cross.mat, par.mat = parent.mat,
                   par.clu = par.clu, Q.eff = Q.eff)

  names(Q.list) <- paste0("Q", 1:length(Q.list))


  ### 2.6 Set the constraint on the different elements and prepare the
  # different elements of the model.


  if ((Q.eff == "cr") || (Q.eff == "biall")){

    trait <- mppData$trait[, 1]

    # 2.6.1 parental or ancestral model

  } else {

    # situations with sum to 0

    if((VCOV == "h.err") && (const == "sum.0")) {

      n.QTL <- length(Q.list)
      trait <- c(mppData$trait[, 1], rep(0, n.QTL))
      cross.mat <- rbind(cross.mat, matrix(0, nrow = n.QTL, ncol = mppData$n.cr))

      add.const <- function(x, Q.mat) {

        Q.mat <- Q.mat[[x]]
        constraint <- matrix(0, nrow = n.QTL, ncol = dim(Q.mat)[2])
        constraint[x, ] <- 1
        Q.mat <- rbind(Q.mat, constraint)
        rownames(Q.mat) <- as.character(1:dim(Q.mat)[1])
        Q.mat

      }

      Q.list <- lapply(X = 1:n.QTL, FUN = add.const, Q.mat = Q.list)

      names(Q.list) <- paste0("Q", 1:length(Q.list))

    # all other VCOV except "h.err" -> set.0 constraint

     } else {

       n.QTL <- length(Q.list)
       trait <- mppData$trait[, 1]

       remove.par <- function(x, QTL, par.ref, par.clu, Q.list){
         
         Q.mat <- Q.list[[x]]
         
         if(Q.eff == "par"){ # parental model
           
           Q.mat[, -which(colnames(Q.mat) == par.ref), drop = FALSE]
           
         } else if (Q.eff == "anc") { # ancestral model
           
           clu.info <- par.clu[which(rownames(par.clu) == QTL[x, 1]), ]
           clu.to.remove <- clu.info[names(clu.info) == par.ref]
           alleles <- as.numeric(substr(colnames(Q.mat), 9,
                                        nchar(colnames(Q.mat))))
           
           Q.mat[, -which(alleles == clu.to.remove), drop = FALSE]
           
         }
         
       }
       
       Q.list <- lapply(X = 1:n.QTL, FUN = remove.par, QTL = QTL,
                        Q.list = Q.list, par.ref = par.ref, par.clu = par.clu)

       names(Q.list) <- paste0("Q", 1:length(Q.list))

   }


 }


  # 3. model computation
  ######################


  model <- QTLModelQeff(mppData = mppData, trait = trait, cross.mat = cross.mat,
                        Q.list = Q.list, VCOV = VCOV)

  # 4. data processing
  ####################

  ### 4.1 sort row model results

  results <- Qeff_res_processing(model = model, mppData = mppData,
                                 cross.mat =  cross.mat, Q.list = Q.list,
                                 QTL = QTL, Q.eff = Q.eff, par.clu = par.clu,
                                 VCOV = VCOV, par.ref = par.ref, const = const)

  ### 4.2 add stars

  stars <- sapply(results[, 4], FUN = sign.star)

  results <- data.frame(results, stars, stringsAsFactors = FALSE)

  colnames(results) <- c("Effect", "Std. Err.", "Test stat.", "p-value", "sign.")
  
  if(VCOV == "h.err"){colnames(results)[3] <- "t-test"} else {
    colnames(results)[3] <- "W-stat"}

  ### 4.3 Complete with parent genotype code

  if(!is.null(mppData$geno.par)){

  results <- add_parents_score(results = results, mppData = mppData,
                               Q.eff = Q.eff, QTL = QTL, par.ref = par.ref,
                               VCOV = VCOV)  

  } 

  ### 4.3 Split into a list with one data.frame for each QTL
  
  
  # Only case where only one term per QTL
  
  if ((Q.eff == "biall") && (is.null(mppData$geno.par))){
    
    partition <- factor(paste0("Q", 1:length(Q.list)),
                        levels = paste0("Q", 1:length(Q.list)))
    
    results <- split(x = results, f = partition)
    
  } else {
    
    # division per cross or per parents.
    
    if(Q.eff == "cr"){
      
      partition <- factor(rep(paste0("Q", 1:length(Q.list)),
                              each = mppData$n.cr),
                          levels = paste0("Q", 1:length(Q.list)))
      
      results <- split(x = results, f = partition)
      
    } else {
      
      partition <- factor(rep(paste0("Q", 1:length(Q.list)),
                              each = mppData$n.par),
                          levels = paste0("Q", 1:length(Q.list)))
      
      results <- split(x = results, f = partition)
      
    }
    
  }
  
  
  return(results)
  
  
}
