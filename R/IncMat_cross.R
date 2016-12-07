################
# IncMat_cross #
################

#' Cross effect incidence matrix
#' 
#' Forms a cross effect incidence matrix.
#' 
#' 
#' @param cross.ind \code{Character} or \code{factor} vector specifying to
#' which cross the genotypes belong.
#' 
#' @return Return:
#' 
#' \item{cross.mat}{Cross effect incidence matrix composed of 0 and 1 value
#' indicating to which cross each genotype belong}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{IncMat_parent}}, \code{\link{IncMat_QTL}}
#' 
#' @examples
#'
#' data(USNAM_mppData)
#' 
#' cross.mat <- IncMat_cross(cross.ind = USNAM_mppData$cross.ind)
#'
#' @export
#'


IncMat_cross <- function(cross.ind){
  
  Cr <- factor(x = cross.ind, levels = unique(cross.ind))
  cross.mat <- model.matrix( ~ Cr - 1, data = Cr)
  
  return(cross.mat)
  
}