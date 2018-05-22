####################
# mppData_chgPheno #
####################

#' Change phenotypic values in a mppData object
#'
#' Replaces the phenotypic values of a \code{mppData} object by new trait values.
#' 
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#' 
#' @param trait Two columns \code{data.frame} with : 1) \code{character}
#' genotypes identifiers; 2) \code{numeric} trait values. \strong{The genotypes
#' identifiers must be identical to \code{mppData$geno.id}.}
#' 
#' @return Return:
#' 
#' \item{mppData}{New \code{mppData} object with replaced phenotypic values.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mppData_subset}},
#' \code{\link{pedigree_update.mppData}}
#' 
#' @examples
#' 
#' # This function is deprecated
#' 
#' # data(USNAM_mppData)
#' # data(USNAM_pheno)
#' 
#' # trait <- USNAM_pheno
#' 
#' # mppData <- mppData_chgPheno(mppData = USNAM_mppData, trait = trait)
#' 

mppData_chgPheno <- function(mppData, trait){
  
  if(!is_mppData(mppData)) {
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  if(!identical(trait[, 1], mppData$geno.id )) {
    
    stop("the genotypes identifiers of 'mppData' and 'trait' are not identical")
    
  }
  
  # substitute the trait
  
  mppData$trait <- subset(x = trait, select = 2, drop = FALSE)
  
  return(mppData)
  
}