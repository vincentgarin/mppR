##############
# IncMat_QTL #
##############

#' QTL incidence matrix
#'
#' Builds a single position QTL incidences matrix
#' 
#' @param x \code{Integer} value indicating the genetic position on the map
#' \code{mppData$map} for which the QTL incidence matrix should be returned.
#' 
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. for more details see \code{\link{mpp_SIM}}. Default = "cr".
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
#' @param cross.mat Cross effect incidence matrix that can be obtained with
#' \code{\link{IncMat_cross}}.
#' 
#' @param par.mat Parents incidence matrix that can be obtained with
#' \code{\link{IncMat_parent}}.
#' 
#' 
#' @return Return:
#' 
#' \item{QTL.mat}{QTL incidence matrix or vector (for bi-allelic model). For the
#' cross-specific model, it represent the difference between the
#' number of allele from parent 2 or B and parent 1 or A divided by two. For
#' parental (ancestral) model it represent the expected number of parent
#' (ancestor) allele copies. For the bi-allelic model, it represent the number
#' of copies of the allele with the highest frequency.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{IncMat_cross}}, \code{\link{IncMat_parent}}, 
#' \code{\link{mpp_SIM}}, \code{\link{parent_cluster}},
#' \code{\link{USNAM_parClu}}
#'
#' @examples
#' 
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' data(USNAM_parClu)
#' 
#' cross.mat <- IncMat_cross(cross.ind = USNAM_mppData$cross.ind)
#' par.mat <- IncMat_parent(USNAM_mppData)
#' par.clu <- USNAM_parClu
#' 
#' QTLmatCr <- IncMat_QTL(x = 2, mppData = USNAM_mppData,
#'                        cross.mat = cross.mat, Q.eff = "cr")
#' 
#' QTLmatPar <- IncMat_QTL(x = 2, mppData = USNAM_mppData, par.mat = par.mat,
#'                         Q.eff = "par")
#' 
#' QTLmatAnc <- IncMat_QTL(x = 2, mppData = USNAM_mppData, par.mat = par.mat,
#'                         par.clu = par.clu, Q.eff = "anc")
#' 
#' QTLmatBi <- IncMat_QTL(x = 2, mppData = USNAM_mppData_bi, Q.eff = "biall")
#' 
#' @export
#' 


IncMat_QTL <- function(x, mppData, Q.eff, par.clu, cross.mat, par.mat) {
  
  pos <- unlist(mppData$map[x, c(2,3)])
  
  if(Q.eff == "cr"){
    
    alpha.pred <- mppData$geno$geno[[pos[1]]]$prob[, pos[2], mppData$n.zigo] -
      mppData$geno$geno[[pos[1]]]$prob[, pos[2], 1]
    
    alpha.pred <- t(rep(1, dim(cross.mat)[2])) %x% alpha.pred
    QTL.mat <- cross.mat * alpha.pred
      
      return(QTL.mat)
    
  } else if(Q.eff == "par") {
    
    # get number of alleles from parent A and B
    
    if (mppData$n.zigo == 3){
      
      alleleA <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 1]) +
        mppData$geno$geno[[pos[1]]]$prob[, pos[2], 2]
      
      alleleB <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 3]) +
        mppData$geno$geno[[pos[1]]]$prob[, pos[2], 2]
      
    } else if (mppData$n.zigo == 2){
      
      alleleA <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 1])
      alleleB <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 2])
      
    }
      
      # form parental QTL matrix
      
      PA_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleA) * par.mat$PA
      PB_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleB) * par.mat$PB
      
      QTL.mat <- PA_pos + PB_pos
    
    return(QTL.mat)
    
  } else if (Q.eff == "anc") {
    
    # form a parental matrix. Same as for Q.eff == "par"
    
    if (mppData$n.zigo == 3){
      
      alleleA <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 1]) +
        mppData$geno$geno[[pos[1]]]$prob[, pos[2], 2]
      
      alleleB <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 3]) +
        mppData$geno$geno[[pos[1]]]$prob[, pos[2], 2]
      
    } else if (mppData$n.zigo == 2){
      
      alleleA <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 1])
      alleleB <- (2*mppData$geno$geno[[pos[1]]]$prob[, pos[2], 2])
      
    }
    
    PA_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleA) * par.mat$PA
    PB_pos <- (t(rep(1, dim(par.mat$PA)[2])) %x% alleleB) * par.mat$PB
    
    QTL.mat <- PA_pos + PB_pos
    
    # modify parental matrix according to ancestral matrix (par. clustering)
    
    A.allele <- as.factor(par.clu[x, ])
                         
    A <- model.matrix(~ A.allele - 1)
    
    QTL.mat <- QTL.mat %*% A
    
    return(QTL.mat)
    
    
  } else if (Q.eff == "biall"){
    
    QTL.mat <- subset(x = mppData$geno, select = x, drop = FALSE)
    
    return(QTL.mat)
    
  }
  
}