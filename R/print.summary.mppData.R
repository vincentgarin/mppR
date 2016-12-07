#########################
# print.summary.mppData #
#########################

#' Print summary.mppData object
#' 
#' @param x object of class \code{summary.mppData}
#' 
#' @param ... further arguments passed to or from other methods. 
#' 
#' @examples 
#' 
#' data(USNAM_mppData)
#' sum.mppData <- summary(USNAM_mppData)
#' print(sum.mppData)
#' 
#' @export
#' 

print.summary.mppData <- function(x, ...){

  if(x$biall) {type.pop <- "object of class 'mppData' for bi-allelic models"
  } else {type.pop <- "object of class 'mppData' for cross-specific, parental or ancestral models"}
  
  cat(type.pop, "\n \n")
  cat("\t Type of population: ", x$typePop, "\n \n")
  cat("\t No. Genotypes: ", x$Ngeno, "\n")
  print(x$par.per.cross, quote = FALSE)
  cat("\n")
  cat("\t Phenotype: ", x$phenoName, "\n")
  cat("\t Percent phenotyped: ", x$phenoPer, "\n \n")
  cat("\t Total marker: ", x$mkNb, "\n")
  cat("\t No. markers:", x$mkChr)
  
}
