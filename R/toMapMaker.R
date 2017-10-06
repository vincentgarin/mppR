##############
# toMapMaker #
##############

#' Map into MapMaker format
#' 
#' Transform a map file with three columns (marker identifiers, chromosomes and
#' positions in cM) into MapMaker format.
#' 
#' @param map Three columns \code{data.frame} with: 1) \code{character}
#' marker identifiers; 2) \code{numeric} integer chromosome indicators
#' (1, 2, 3,...); and 3) \code{numeric} positions in centi-Morgan.
#' 
#' @return Return:
#' 
#' \item{new.map}{Map with the MapMaker format: \code{list} of matrix with one
#' element for each chromosome. Each single chromosome matrix is composed of four
#' columns: 1) marker identifiers; 2) dist.to: one lagged difference
#' positions in cM with Inf at first position; 3) dist.from: one lagged
#' difference positions in cM with Inf at last position ; 4) locus: marker
#' positions in cM.}
#' 
#' @author Vincent Garin
#' 
#' @examples
#' 
#' data(USNAM_map)
#' 
#' map_MapMaker <- toMapMaker(USNAM_map)
#' 
#' @export
#' 


toMapMaker <- function(map) {
  
  chr.id <- levels(as.factor(map[, 2])) # chromosome id
  
  chr_list <- vector(mode = "list", length = length(chr.id)) # empty list
  names <- c()
  
  for (i in seq_along(chr.id)) {
    
    # select one chromosome
    
    map.sel <- map[which(map[, 2] == i), ]
    
    distance <- diff(x = map.sel[, 3], lag = 1)
    dist.to <- c(Inf, distance)
    dist.from <- c(distance, Inf)
    chr.mat <- data.frame(map.sel[, 1], dist.to, dist.from, map.sel[, 3],
                          stringsAsFactors = FALSE)
    colnames(chr.mat) <- c("markers", "dist.to", "dist.from", "locus")
    names <- c(names, paste("chr", i, sep = ""))
    
    # fill the list
    
    chr_list[[i]] <- chr.mat
    
  }
  
  names(chr_list) <- names
  
  return(chr_list)
  
}