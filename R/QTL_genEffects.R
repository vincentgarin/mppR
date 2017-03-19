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
#' is possible to estimate maximum n-1 parental or ancestral alleles per connected
#' part of the design. Connected parts of the design can be determined using
#' Weeks and Williams (1964) method (\code{\link{design_connectedness}}).
#' Connected parts of the design are indicated in the results and parental
#' (ancestral) allele score getting 0 value are the one that could not be
#' estimated due to singularities. Estimated effects of the other alleles should
#' be interpreted within connected part as deviation with respect to the
#' non-estimable effect(s).
#' 
#' The user can specify if the allele effect that should be put as reference
#' within connected parts are the most or less used one setting argument
#' \code{ref.all.most = TRUE or FALSE}. Allele usage is first defined
#' in term of number of crosses where the allele segregates. Then, if two
#' alleles segregate in an equal number of crosses, we look at the total cross
#' sizes. For example in a NAM design the central parent is the most used allele.
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
#' @param ref.all.most \code{Logical} value specifying which parental
#' (ancestral) allele should be used as reference. If
#' \code{ref.all.most = TRUE}, within each connected part of the design,
#' the most used allele will be used as reference. Allele usage is first defined
#' in term of number of crosses where the allele segregates. Then, if two
#' alleles segregate in an equal number of crosses, we look at the total cross
#' sizes. If \code{ref.all.most = FALSE}, the less used allele will be used as
#' reference within each connected part. Default = TRUE.
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
#' \item{For parental and ancestral model, indicator of connected part of the
#' design and reference.}
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
                           VCOV = "h.err",  ref.all.most = TRUE) {
  
  # 1. Check data format
  ######################
  
  check.model.comp(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
                   par.clu = par.clu, QTL = QTL, fct = "QTLeffects")
  
  
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
  
  if ((Q.eff == "par") || (Q.eff == "anc")){
    
    if(Q.eff == "par"){
      
      # 1. determine the connected parts
      
      con.part <- design_connectedness(par.per.cross = mppData$par.per.cross,
                                        plot.des = FALSE)
      
      # 2. determine the most (less) used allele within connected part
      
      allele_order <- c()
      allele_ref <- c()
      
      for(i in seq_along(con.part)){
        
        con.part_i <- con.part[[i]]
        
        # subset the par.per.cross object
        
        index <- apply(X = mppData$par.per.cross[, c(2,3)], MARGIN = 1,
                       FUN = function(x, ref) sum(x %in% ref) > 0,
                       ref = con.part_i)
        
        par.per.cross_i <- mppData$par.per.cross[index, , drop = FALSE]
        
        # susbet the cross indicator according to the cross retained
        
        cross.ind_i <- mppData$cross.ind[mppData$cross.ind %in% par.per.cross_i[, 1]]
        
        allele_ord_i <- most.used.allele(par.per.cross_i, cross.ind_i,
                                         most.used = !ref.all.most)
        
        allele_order <- c(allele_order, allele_ord_i)
        allele_ref <- c(allele_ref, allele_ord_i[length(allele_ord_i)])
        
      }
      
      # order the parental alleles in QTL incidence matrices
      
      Q.list <- lapply(X = Q.list, FUN = function(x, ref) x[, ref],
                       ref = allele_order)
      
      
    } else if (Q.eff == "anc"){
      
      if(is.character(QTL)){ QTL.list <- QTL} else {QTL.list <- QTL[, 1]}
      
      con.part <- vector(mode = "list", length(Q.list))
      allele_ref <- vector(mode = "list", length(Q.list))
      
      for(i in 1:length(Q.list)){
        
        # 1. Determine the connected parts
        
        par.clu_i <- par.clu[QTL.list[i], ]
        par.clu_i <- paste0("A.allele", par.clu_i)
        names(par.clu_i) <- mppData$parents
        
        all.p1 <- par.clu_i[mppData$par.per.cross[, 2]]
        all.p2 <- par.clu_i[mppData$par.per.cross[, 3]]
        
        par.per.cross_i <- cbind(mppData$par.per.cross[, 1], all.p1, all.p2)
        
        con.part_i <- design_connectedness(par.per.cross = par.per.cross_i,
                                            plot.des = FALSE)
        
        con.part[[i]] <- con.part_i
        
        # 2. Within the connected parts determine the most (least) used allele
        
        allele_ord_i <- c()
        allele_ref_i <- c()
        
        for(j in seq_along(con.part_i)){
          
          con.part_j <- con.part_i[[j]]
          
          # subset the par.per.cross object
          
          index <- apply(X = par.per.cross_i[, c(2, 3)], MARGIN = 1,
                         FUN = function(x, ref) sum(x %in% ref) > 0,
                         ref = con.part_j)
          
          par.per.cross_j <- par.per.cross_i[index, , drop = FALSE]
          
          # susbet the cross indicator according to the cross retained
          
          cross.ind_j <- mppData$cross.ind[mppData$cross.ind %in% par.per.cross_j[, 1]]
          
          allele_ord_j <- most.used.allele(par.per.cross_j, cross.ind_j,
                                           most.used = !ref.all.most)
          
          allele_ord_i <- c(allele_ord_i, allele_ord_j)
          allele_ref_i <- c(allele_ref_i, allele_ord_j[length(allele_ord_j)])
          
        }
        
        # Store the results
        
        allele_ref[[i]] <- allele_ref_i
        
        # 3. Order the QTL incidence matrix
        
        Q.list[[i]] <- Q.list[[i]][, allele_ord_i]
        
      }
      
    }
    
    
  } else {
    
    con.part <- allele_ref <- NULL
    
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
                                 VCOV = VCOV, allele_ref = allele_ref,
                                 con.part = con.part)
  
  return(results)
  
}