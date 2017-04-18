##########
# MQE_R2 #
##########

#' Global and partial R squared for multi-QTL effects 
#' 
#' Computes the global and partial (adjusted) R squared of a list of QTLs
#' including positions with different type of QTL effects. For example
#' position one is a parental effect, position two a bi-allelic QTL, etc.
#' The type of effect of the QTL position are specified in \code{Q.eff}.
#' 
#' The R squared computation is done using a linear model
#' corresponding to \code{VCOV = 'h.err'} For more details about R squared
#' computation and adjustement look at function \code{\link{QTL_R2}}.
#'
#' @param mppData An IBD object of class \code{mppData}
#' See \code{\link{mppData_form}} for details. Default = NULL.
#'
#' @param mppData_bi Required IBS object of class \code{mppData} if the user
#' wants to allow QTLs with a bi-allelic effect. \strong{The list of marker must
#' be strictly the same as the one of \code{mppData}.} Default = NULL.
#' 
#' @param QTL Vector of \code{character} markers or inbetween marker positions
#' names. Default = NULL.
#' 
#' @param Q.eff \code{Character} vector indicating for each QTL position the
#' type of QTL effect among: "cr", "par", "anc" and "biall". For details look at
#' \code{\link{mpp_SIM}}.
#'
#' @param par.clu Required argument if the user wants to allow QTLs with an
#' ancestral effect. \code{interger matrix} representing the results of a parents
#' genotypes
#' clustering. The columns represent the parental lines and the rows
#' the different markers or in between positions. \strong{The columns names must
#' be the same as the parents list of the mppData object. The rownames must be
#' the same as the map marker list of the mppData object.} At a particular
#' position, parents with the same value are assumed to inherit from the same
#' ancestor. for more details, see \code{\link{USNAM_parClu}} and
#' \code{\link{parent_cluster}}. Default = NULL.
#'
#' @param glb.only \code{Logical} value. If \code{LR.R2 = TRUE}, only the global
#' and global adjusted R squared will be returned. Default = FALSE.
#' 
#'
#' @return Return:
#'
#' List containing the following objects:
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
#' @seealso \code{\link{mppData_form}}, \code{\link{mpp_SIM}},
#' \code{\link{parent_cluster}},
#' \code{\link{QTL_R2}}, \code{\link{USNAM_parClu}}
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
#' QTL <- QTL_select(SIM)
#' QTL <- QTL[, 1]
#' 
#' MQE_R2(mppData = mppData, mppData_bi = mppData_bi, QTL = QTL,
#'        Q.eff = c("par", "anc", "biall"), par.clu = par.clu)
#'
#' @export
#'


MQE_R2 <- function(mppData = NULL, mppData_bi = NULL, QTL = NULL, Q.eff,
                   par.clu = NULL, glb.only = FALSE){
  
  # 1. check the data format
  ##########################
  
  check.MQE(mppData = mppData, mppData_bi =  mppData_bi, Q.eff = Q.eff,
            par.clu = par.clu, QTL = QTL, fct = "R2")
  
  # 2. elements for the model
  ###########################
  
  ### 2.1 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.2 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.3 modify the par.clu object order parents columns and replace
  # monomorphic
  
  if ("anc" %in% Q.eff) {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.4 Formation of the list of QTL incidence matrices
  
  # order list of QTL positions
  
  Q.pos <- vapply(X = QTL,
                  FUN = function(x, mppData) which(mppData$map[, 1] == x),
                  FUN.VALUE = numeric(1), mppData = mppData)
  
  Q.ord <- data.frame(QTL, Q.eff, Q.pos, stringsAsFactors = FALSE)
  
  Q.ord <- Q.ord[order(Q.pos), ]
  
  QTL <- Q.ord[, 1]; Q.eff <- Q.ord[, 2]
  
  Q.pos <- which(mppData$map[, 1] %in% QTL)
  
  # form a list of QTL incidence matrices with different type of QTL effect.
  
  # function to produce different type of QTL incidence matricdes
  
  
  Q.list <- mapply(FUN = IncMat_QTL_MQE, x = Q.pos, Q.eff = Q.eff,
                   MoreArgs = list(mppData = mppData, mppData_bi = mppData_bi,
                                   par.clu = par.clu, cross.mat = cross.mat,
                                   par.mat = parent.mat, order.MAF = TRUE),
                   SIMPLIFY = FALSE)
  
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
      
      ### 3.2 Compute the partial R squared
      
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