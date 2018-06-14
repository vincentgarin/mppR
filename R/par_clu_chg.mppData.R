#######################
# par_clu_chg.mppData #
#######################

#' Change the parent clustering in a \code{mppData} object
#'
#' Change the parent clustering in a \code{mppData} object.
#' 
#' @param mppData An object of class \code{mppData}. The \code{mppData} object
#' must have been processed with \code{\link{parent_cluster.mppData}}
#' 
#' @param par.clu \code{Interger matrix} representing the results of a
#' parents genotypes clustering. The columns represent the parental lines and
#' the rows the markers. The columns names must be the same as the parents
#' list of the mppData object. The rownames must be the same as the map marker
#' list of the mppData object. At a particular position, parents with the same
#' value are assumed to inherit from the same ancestor. for more details,
#' see \code{\link{par_clu}}.
#' 
#' @return 
#' 
#' \code{mppData} object with the new parent clustering information.
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{parent_cluster.mppData}}, \code{\link{par_clu}}
#' 
#' @examples
#' 
#' data(mppData)
#' data(par_clu)
#' 
#' mppData <- par_clu_chg.mppData(mppData = mppData, par.clu = par_clu)
#' 
#' 
#' @export
#' 

par_clu_chg.mppData <- function(mppData, par.clu){
  
  # 1. Check argument format
  ##########################
  
  ### mppData
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # test if correct step in the mppData processing
  
  if(mppData$status != 'complete'){
    
    stop("'mppData' must have been processed with parent_cluster.mppData")
    
  }
  
  ### par.clu
  
  if(!is.matrix(par.clu)){
    
    stop("'par.clu' argument is not a matrix")
    
  }
  
  if(!is.integer(par.clu)){
    
    stop("'par.clu' is not integer")
    
  }
  
  # list parent
  
  new_par <- colnames(par.clu)
  
  if(!all(new_par %in% mppData$parents)) {
    
    wrong.par <- new_par[!(new_par %in% mppData$parents)]
    pbpar <- paste(wrong.par, collapse = ", ")
    
    message <- sprintf(ngettext(length(wrong.par),
                                "the following parent %s is not present in the old 'par.clu'",
                                "the following parents %s are not present in the old 'par.clu'"),
                       pbpar)
    
    stop(message)
    
  }
  
  # list markers
  
  if(!identical(rownames(par.clu), mppData$map[, 1])){
    
    stop("the markers of 'par.clu' and 'mppData' are not identical")
    
  }
  
  # 2. check monomorphism in par.clu
  ##################################
  
  par.clu <- par.clu[, mppData$parents]
  
  par_clu <- parent_clusterCheck(par.clu = par.clu)
  
  # 3. Recalculate the number of ancestral cluster
  ################################################
  
  nb.cl <- apply(X = par_clu[[1]], MARGIN = 1,
                 FUN = function(x) length(unique(x)))
  
  av.cl <- mean(nb.cl)
  
  mppData$par.clu <- par_clu[[1]]
  
  mppData$n.anc <- av.cl
  
  return(mppData)
  
}