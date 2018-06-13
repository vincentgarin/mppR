###########################
# pedigree_update.mppData #
###########################

#' Modify pedigree information in a \code{mppData} object
#'
#' Modify the pedigree information of a \code{mppData} object.
#' 
#' In the formation of the \code{mppData} object, the functions only
#' uses the parental relationships specified in the \code{par.per.cross}
#' argument. However, if pedigree relationship among parents are know by the
#' user a more complete pedigree information can be provided to the
#' \code{mppData} object. This pedigree information can be given via the
#' \code{pedigree} argument of the \code{pedigree_update.mppData()} function. 
#' 
#' @param mppData An object of class \code{mppData}. The \code{mppData} object
#' must have been processed with \code{\link{QC.mppData}}
#' 
#' @param pedigree Four columns \code{data.frame}: 1) the type of
#' genotype: "offspring" for the last genration and "founder" for genotype above
#' the offspring in the pedigree; 2) the genotype indicator of offsrping and
#' founders; 3-4) the parent 1 (2) of each line. \strong{The genotypes
#' identifiers of the offspring must be identical to \code{mppData$geno.id}.
#' The row giving the pedigree of an individual must appear before any row
#' where that individual appears as a parent.}
#' 
#' @return
#' 
#' \code{mppData} object with completed pedigree information.
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{QC.mppData}}, \code{\link{subset.mppData}}
#' 
#' @examples
#' 
#' data(mppData)
#' pedigree <- mppData$ped.mat
#' 
#' # Let us assume for example that parent CML103 and CML322 share a
#' # common parent (A1). Then we can write a new pedigree adding this
#' # information.
#' 
#' founder.info <- data.frame(rep("founder", 2), c("CML103", "CML322"),
#'                            rep("A1", 2), c("A2", "A3"))
#' colnames(founder.info) <- colnames(pedigree)
#' 
#' pedigree.new <- rbind(founder.info, pedigree)
#' 
#' mppData <- pedigree_update.mppData(mppData = mppData,
#' pedigree = pedigree.new)
#' 
#' 
#' @export
#' 


pedigree_update.mppData <- function(mppData, pedigree){
  
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  if(mppData$status == 'init'){
    
    stop("'mppData' must at least have been processed with QC.mppData")
    
  }
  
  if(!identical(as.character(pedigree[pedigree[, 1] == "offspring", 2]),
                mppData$geno.id )) {
    
    stop("the offspring genotypes identifiers of 'mppData' and 'pedigree' are not identical")
    
  }
  
  # Test if the newly introduced genotype are also mentioned earlier
  
  new.geno <- pedigree[pedigree[, 1] == "founder", 2]
  
  old.geno <- unlist(pedigree[pedigree[, 1] == "offspring", 3:4])
  
  if(!all(new.geno %in% old.geno)) {
    
    wrong.geno <- new.geno[!(new.geno %in% old.geno)]
    pbgeno <- paste(wrong.geno, collapse = ", ")
    
    message <- sprintf(ngettext(length(wrong.geno),
                                "'pedigree' genotype %s is not present in 'mppData'",
                                "'pedigree' genotypes %s are not present in 'mppData'"),
                       pbgeno)
    
    stop(message)
    
  }
  
  mppData$ped.mat <- pedigree
  
  return(mppData)
  
}