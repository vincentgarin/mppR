##################
# create.mppData #
##################

#' Create a multi-parent population data object
#' 
#' This function combines all raw data sources in a single data object of class
#' \code{mppData}.
#' 
#' @param geno.off \code{Character} marker score \code{matrix} of the offspring
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the offspring genotypes identifiers similar to
#' the one used in \code{pheno}. The columns names must be the marker names
#' similar to the one used in \code{map}. Marker scores must be coded using one
#' letter per allele. For example, AA, CC, GG, etc. Missing values must be coded
#' \code{NA}.}
#' 
#' @param geno.par \code{Character} marker score \code{matrix} of the parents
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the parents genotypes identifiers similar to
#' the one used in \code{par.per.cross}. The columns names must be the marker
#' names similar to the one used in \code{map}. Marker scores must be coded
#' using one letter per allele. For example, AA, CC, GG, etc. Missing values
#' must be coded \code{NA}.}
#' 
#' @param map Three columns \code{data.frame} with: 1) marker or identifiers;
#' 2) chromosome; 3) positions in centi-Morgan.\strong{
#' The marker identifiers must be identical to the column names of the maker
#' matrices (\code{geno.off} and \code{geno.par}).}
#' 
#' @param pheno \code{Numeric matrix} with one column per trait and rownames
#' as genotpes identifiers. \strong{The genotypes identifiers must be identical
#' to the rownames of the offspring marker matrix (\code{geno.off}).}
#' 
#' @param cross.ind \code{Character} vector indicating to which cross each
#' genotype belongs.
#' 
#' @param par.per.cross Three columns \code{Character matrix} specifying :
#' 1) the cross indicators; 2) the parents 1 identifiers
#' of the crosses; 3) the parents 2 identifiers of the crosses.
#' \strong{The list of crosses must contain the same cross indicators as in
#' \code{cross.ind} and they must appear in the same order.
#' The list of parent identifiers must be the same to the rownames of
#' the argument \code{geno.par}}.
#' 
#' @return 
#' 
#' a \code{list} of class \code{mppData} which contains the following elements
#' 
#' \item{geno.off}{\code{Matrix} with offspring marker scores.}
#' 
#' \item{geno.par}{\code{Matrix} with parents marker scores.}
#' 
#' \item{pheno}{\code{Matrix} with phenotypic trait values.}
#' 
#' \item{map}{\code{Data.frame} with genetic marker information.}
#' 
#' \item{cross.ind}{Cross indicator.}
#' 
#' \item{par.per.cross}{\code{Character matrix} information about cross and
#' the parents of the crosses.}
#' 
#' The \code{list} also contain other arguments that will be filled later in
#' the data processing.
#' 
#' @author Vincent Garin
#' 
#' @examples
#' 
#' data(USNAM_geno)
#' data(USNAM_map)
#' data(USNAM_pheno)
#' 
#' geno.off <- USNAM_geno[7:506, ]
#' geno.par <- USNAM_geno[1:6, ]
#' map <- USNAM_map
#' pheno <-  USNAM_pheno
#' cross.ind <- substr(rownames(pheno), 1, 4)
#' par.per.cross <- cbind(unique(cross.ind), rep("B73", 5),
#'                        rownames(geno.par)[2:6])
#' 
#' mppData <- create.mppData(geno.off = geno.off, geno.par = geno.par,
#'                           map = map, pheno = pheno, cross.ind = cross.ind,
#'                           par.per.cross = par.per.cross)
#' 
#' @export
#' 

create.mppData <- function(geno.off = NULL, geno.par = NULL, map = NULL,
                           pheno = NULL, cross.ind = NULL,
                           par.per.cross = NULL){
  
  # 1. Check the correct format of the raw data
  #############################################
  
  # Test if one object is missing
  
  if(is.null(geno.off)){stop('Argument geno.off is not provided.')}
  
  if(is.null(geno.par)){stop('Argument geno.par is not provided.')}
  
  if(is.null(map)){stop('Argument map is not provided.')}
  
  if(is.null(pheno)){stop('Argument pheno is not provided.')}
  
  if(is.null(cross.ind)){stop('Argument cross.ind is not provided.')}
  
  if(is.null(par.per.cross)){stop('Argument par.per.cross is not provided.')}
  
  # test format of the marker matrices (geno.off, geno.par)
  
  if(!is.matrix(geno.off)){ stop("The geno.off argument is not a matrix.") }
  
  if(!is.character(geno.off)){
    
    stop("The geno.off argument is not a character matrix.")
    
  }
  
  if(!is.matrix(geno.par)){ stop("The geno.par argument is not a matrix.") }
  
  if(!is.character(geno.par)){
    
    stop("The geno.par argument is not a character matrix.")
    
  }
  
  # test if the marker identifiers are the same in the map and in the marker
  # matrices
  
  if(!identical(colnames(geno.off), map[, 1])){
    
    stop(paste("The marker identifier of the offspring marker matrix (geno.off)",
               "and the list of marker (map[, 1]) are not the same."))
    
  }
  
  if(!identical(colnames(geno.par), map[, 1])){
    
    stop(paste("The marker identifier of the parent marker matrix (geno.par)",
               "and the list of marker (map[, 1]) are not the same."))
    
  }
  
  # test phenotype values
  
  if(!is.matrix(pheno)){stop('The pheno argument is not a matrix.')}
  
  if(!is.numeric(pheno)){stop('The pheno argument is not a numeric matrix.')}
  
  # test if the list of genotype is the same between the offspring marker matrix
  # and the phenotypic values.
  
  if(!identical(rownames(geno.off), rownames(pheno))){
    
    stop(paste("The genotype identifiers of the offspring marker matrix",
               "(geno.off) and the list of genotypes (rownames(pheno)) are not",
               "the same."))
    
  }
  
  # length cross.ind same as list of genotypes
  
  if(length(cross.ind) != dim(geno.off)[1]){
    
    stop(paste("The cross indicator vector (cross.ind) length does not have",
               "the same length as the genotype list (dim(geno.off)[1])"))
  }
  
  # test par.per.cross
  
  if(!is.matrix(par.per.cross)){
    
    stop("The par.per.cross argument is not a matrix.")
    
  }
  
  if(!is.character(par.per.cross)){
    
    stop("The par.per.cross argument is not a character matrix.")
    
  }
  
  # remove the eventual rownames of par.per.cross
  
  if(!is.null(rownames(par.per.cross))){
    
    rownames(par.per.cross) <- NULL
    
  }
  
  if (!identical(unique(cross.ind), par.per.cross[, 1])){
    
    stop(paste("The cross identifiers used in cross.ind and",
               "in par.per.cross are different."))
    
  }
  
  # test the similarity of parents list between par.per.cross and
  # rownames(geno.par)
  
  parents <- union(par.per.cross[, 2], par.per.cross[, 3])
  
  if(sum(!(parents %in% rownames(geno.par)))>0){
    
    list.par <- paste(parents[!(parents %in% rownames(geno.par))])
    
    stop(paste("The following parents indicators:", list.par,
               "(is) are present in par.per.cross object but not in the",
               "rownames of the geno.par matrix"))
    
  }
  
  # Check the number of connected parts
  
  con_part <- design_connectivity(par.per.cross, plot_des = FALSE)
  n_con <- length(con_part)
  
  n_geno <- dim(geno.off)[1]
  n_par <- dim(geno.par)[1]
  n_cr <- dim(par.per.cross)[1]
  n_pheno <- dim(pheno)[2]
  
  # Possibility to check here the necessity to have at least 2 crosses and 3
  # different parents.
  
  ####### end check format
  
  mppData <- list(geno.off = geno.off, geno.IBS = NULL, geno.IBD = NULL,
                  geno.id = NULL, ped.mat = NULL, allele.ref = NULL,
                  geno.par = geno.par, geno.par.clu = NULL, par.clu = NULL,
                  pheno = pheno, map = map, haplo.map = NULL,
                  cross.ind = cross.ind, par.per.cross = par.per.cross,
                  type = NULL, parents = NULL, n.cr = NULL, n.par = NULL,
                  n.zigo = NULL, rem.mk = NULL, rem.gen = NULL,
                  status = 'init')
  
  class(mppData) <- c("mppData", "list")
  
  # Return a short message that everything was successful
  
  cat("\n")
  cat("mppData object created!\n")
  cat("\n")
  cat(paste(n_geno, 'genotypes\n'))
  cat(paste(n_cr, 'crosses\n'))
  cat(paste(n_par, 'parents\n'))
  cat(paste(n_pheno, 'phenotype(s)\n'))
  if(n_con == 1){
    
    cat('1 connected part\n')
    
  } else {message(paste(n_con, 'independent connected parts\n'))}
  
  
  return(mppData)
  
}