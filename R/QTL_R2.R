##########
# QTL_R2 #
##########

#' QTL global and partial R squared
#' 
#' Computes the global and partial (adjusted) R squared of a  list of QTLs using
#' a linear model.
#' 
#' The function computes R squared statistics using a linear model. This model
#' correspond to a model where \code{VCOV = "h.err"} (for details see
#' \code{\link{mpp_SIM}}). The extra variance explained by a full model
#' containing the QTL terms with respect to a reduced model containing only the
#' cross intercept terms and uses the ratio between the residual sum of square
#' of these two models: R2 = 1-(RSS(f))/(RSS(r)).
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
#' @param glb.only \code{Logical} value. If glb.only = TRUE, only the global and 
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
#' QTL_R2(mppData = USNAM_mppData, QTL = QTL, Q.eff = "cr")
#' 
#' 
#' @export
#'


QTL_R2 <- function(mppData, QTL = NULL, Q.eff = "cr", par.clu = NULL,
                   glb.only = FALSE){
  
  # 1. check the data format
  ##########################
  
  check.model.comp(mppData = mppData, Q.eff = Q.eff, par.clu = par.clu,
                   QTL = QTL, fct = "R2")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.2 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.3 modify the par.clu object order parents columns and replace
  # monomorphic
  
  if (Q.eff == "anc") {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.4 Formation of the list of QTL
  
  if(is.character(QTL)){
    
    Q.pos <- which(mppData$map[, 1] %in% QTL)
    
  } else {
    
    Q.pos <- which(mppData$map[, 1] %in% QTL[, 1])
    
  }
  
  Q.list <- lapply(X = Q.pos, FUN = IncMat_QTL, mppData = mppData,
                   cross.mat = cross.mat, par.mat = parent.mat,
                   par.clu = par.clu, Q.eff = Q.eff, order.MAF = TRUE)
  
  n.QTL <- length(Q.list)
  
  # 3. Compute the R squared
  ##########################
  
  ### 3.1 Global adjusted and unadjusted linear R squared
  
  R2.all <- R2_lin(mppData = mppData, QTL = do.call(cbind, Q.list))
  
  R2 <- R2.all[[1]]
  R2.adj <- R2.all[[2]]
  
  
  if(glb.only) {
    
    return(list(glb.R2 = R2, glb.adj.R2 = R2.adj))
    
  } else {
    
    if(n.QTL > 1){
      
      # functions to compute the R squared or all QTL minus 1 or only 1 QTL position
      
      part.R2.diff <- function(x, QTL, mppData) {
        R2_lin(mppData = mppData, QTL = do.call(cbind, Q.list[-x]))
      }
      
      part.R2.sg <- function(x, QTL, mppData) {
        R2_lin(mppData = mppData, QTL = do.call(cbind, Q.list[x]))
      }
      
      R2.dif <- lapply(X = 1:n.QTL, FUN = part.R2.diff, QTL = Q.list,
                       mppData = mppData)
      
      R2_i.dif <- lapply(X = R2.dif, FUN = function(x) x[[1]])
      R2_i.dif.adj <- lapply(X = R2.dif, FUN = function(x) x[[2]])
      
      R2_i.dif <- R2 - unlist(R2_i.dif) # difference full model and model minus i
      R2_i.dif.adj <- R2.adj - unlist(R2_i.dif.adj)
      
      R2.sg <- lapply(X = 1:n.QTL, FUN = part.R2.sg, QTL = Q.list,
                      mppData = mppData)
      
      R2_i.sg <- unlist(lapply(X = R2.sg, FUN = function(x) x[[1]]))
      R2_i.sg.adj <- unlist(lapply(X = R2.sg, FUN = function(x) x[[2]]))
      
      
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