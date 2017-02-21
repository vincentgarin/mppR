###################
# summary.mppData #
###################

#' Print summary of mppData object
#' 
#' S3 \code{summary} method for oject of class \code{mppData}. The summary of
#' the map information only list the marker and not the added in between
#' positions.
#' 
#' @param object An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#' 
#' @param ... further arguments passed to or from other methods.
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{USNAM_mppData}},
#' \code{\link{USNAM_mppData_bi}}
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' summary(USNAM_mppData)
#' 
#' @export
#' 


summary.mppData <- function(object, ...) {
  
  stopifnot(inherits(object, "mppData"))
  ans <- list()
  
  par.per.cross <- t(object$par.per.cross)
  colnames(par.per.cross) <- attr(table(object$cross.ind), "dimnames")[[1]]
  par.per.cross <- rbind(par.per.cross, table(object$cross.ind))
  colnames(par.per.cross) <- rep(" ", dim(par.per.cross)[2])
  rownames(par.per.cross) <- c("Crosses", "Parent 1", "Parent 2", "N")
  
  # genotype info
  
  ans$typePop <- object$type
  ans$Ngeno <- length(object$geno.id)
  ans$par.per.cross <- par.per.cross
  ans$NgenoCr <- table(object$cross)
  ans$biall <- object$biall
  
  # phenotype info
  
  ans$phenoName <- colnames(object$trait)
  ans$phenoPer <- (1 - (sum(is.na(object$trait[, 1]))/
                          length(object$trait[, 1])))*100
  
  # map information
  
  if(!object$biall){
    
    #remove the inbetween positions of the map
    map <- object$map
    inb.pos.id <- paste0(unique(map$chr), "_loc")
    inb.pos.loc <- data.frame(lapply(X = inb.pos.id,
                                     FUN = function(y) grepl(pattern = y, x = map[, 1])))
    map <- map[rowSums(x = inb.pos.loc) == 0, ]
    
    ans$mkNb <- dim(map)[1]
    ans$mkChr <- table(map[, 2])
    
    
  } else {
    ans$mkNb <- dim(object$map)[1]
    ans$mkChr <- table(object$map[, 2])
  }
  
  class(ans) <- "summary.mppData"
  ans
  
  }