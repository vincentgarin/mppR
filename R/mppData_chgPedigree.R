#######################
# mppData_chgPedigree #
#######################

#' Change pedigree information in a mppData object
#'
#' Replaces the pedigree information of a \code{mppData} object.
#' 
#' In the formation of the \code{mppData} object\code{mppData_from()} only
#' uses the parental relationships specified in the \code{par.per.cross}
#' argument. However, if pedigree relationship among parents are know by the
#' user a more complete pedigree information can be provided to the
#' \code{mppData} object. This pedigree information can be given via the
#' \code{pedigree} argument of the \code{mppData_chgPedigree()} function. 
#' 
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#' 
#' @param pedigree Four columns \code{data.frame}: 1) the type of
#' genotype: "offspring" for the last genration and "founder" for genotype above
#' the offspring in the pedigree; 2) the genotype indicator of offsrping and
#' founders; 3-4) the parent 1 (2) of each line. \strong{The genotypes
#' identifiers of the offspring must be identical to \code{mppData$geno.id}.
#' The row giving the pedigree of an individual must appear before any row
#' where that individual appears as a parent.}
#' 
#' @return Return:
#' 
#' \item{mppData}{New \code{mppData} object with replaced pedigree information}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{mppData_subset}},
#' \code{\link{mppData_chgPheno}}
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' pedigree <- USNAM_mppData$ped.mat
#' 
#' # Let us assume for example that parent CML103 and CML322 share a
#' # common parent (A1). Then we can write a new pedigree adding this
#' # information on the top.
#' 
#' founder.info <- data.frame(rep("founder", 2), c("CML103", "CML322"),
#'                            rep("A1", 2), c("A2", "A3"))
#' colnames(founder.info) <- colnames(pedigree)
#' 
#' pedigree.new <- rbind(founder.info, pedigree)
#' 
#' mppData <- mppData_chgPedigree(mppData = USNAM_mppData,
#' pedigree = pedigree.new)
#' 
#' 
#' @export
#' 


mppData_chgPedigree <- function(mppData, pedigree){
  
  
  stopifnot(inherits(mppData, "mppData"))
  
  if(!identical(as.character(pedigree[pedigree[, 1]=="offspring", 2]),
                mppData$geno.id )) {
    
    stop("The genotypes  identifiers of the mppData object and the offspring
         pedigree information are not the same")
    
  }
  
  mppData$ped.mat <- pedigree
  
  return(mppData)
  
}