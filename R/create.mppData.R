##################
# create.mppData #
##################

#' Create a multi-parent population data object
#' 
#' This function combines all raw data sources in a single data object of class
#' \code{mppData}.
#' 
#' @param geno_off \code{Character} marker score \code{matrix} of the offspring
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the offspring genotypes identifiers similar to
#' the one used in \code{pheno}. The columns names must be the marker names
#' similar to the one used in \code{map}. Marker scores must be coded using one
#' letter per allele. For example, AA, CC, GG, etc. Missing values must be coded
#' \code{NA}.}
#' 
#' @param geno_par \code{Character} marker score \code{matrix} of the parents
#' with genotypes as row and markers as column.
#' \strong{Rows names must be the parents genotypes identifiers similar to
#' the one used in \code{par_per_cross}. The columns names must be the marker
#' names similar to the one used in \code{map}. Marker scores must be coded
#' using one letter per allele. For example, AA, CC, GG, etc. Missing values
#' must be coded \code{NA}.}
#' 
#' @param map Three columns \code{data.frame} with: 1) marker or identifiers;
#' 2) chromosome; 3) positions in centi-Morgan.\strong{
#' The marker identifiers must be identical to the column names of the maker
#' matrices (\code{geno_off} and \code{geno_par}).}
#' 
#' @param pheno \code{Numeric matrix} with one column per trait and rownames
#' as genotpes identifiers. \strong{The genotypes identifiers must be identical
#' to the rownames of the offspring marker matrix (\code{geno_off}).}
#' 
#' @param cross_ind \code{Character} vector indicating to which cross each
#' genotype belongs.
#' 
#' @param par_per_cross Three columns \code{Character matrix} specifying :
#' 1) the cross indicators; 2) the parents 1 identifiers
#' of the crosses; 3) the parents 2 identifiers of the crosses.
#' \strong{The list of crosses must contain the same cross indicators as in
#' \code{cross_ind} and they must appear in the same order.
#' The list of parent identifiers must be the same to the rownames of
#' the argument \code{geno_par}}.
#' 
#' @return 
#' 
#' a \code{list} of class \code{mppData} which contains the following elements
#' 
#' \item{geno_off}{\code{Matrix} with offspring marker scores.}
#' 
#' \item{geno_par}{\code{Matrix} with parents marker scores.}
#' 
#' \item{pheno}{\code{Matrix} with phenotypic trait values.}
#' 
#' \item{map}{\code{Data.frame} with genetic marker information.}
#' 
#' \item{cross_ind}{Cross indicator.}
#' 
#' \item{par_per_cross}{\code{Character matrix} information about cross and
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
#' geno_off <- USNAM_geno[7:506, ]
#' geno_par <- USNAM_geno[1:6, ]
#' map <- USNAM_map
#' pheno <-  USNAM_pheno
#' cross_ind <- substr(rownames(pheno), 1, 4)
#' par_per_cross <- cbind(unique(cross_ind), rep("B73", 5),
#'                        rownames(geno_par)[2:6])
#' 
#' mppData <- create.mppData(geno_off = geno_off, geno_par = geno_par,
#'                           map = map, pheno = pheno, cross_ind = cross_ind,
#'                           par_per_cross = par_per_cross)
#' 
#' @export
#' 

create.mppData <- function(geno_off = NULL, geno_par = NULL, map = NULL,
                           pheno = NULL, cross_ind = NULL,
                           par_per_cross = NULL){
  
  # 1. Check the correct format of the raw data
  #############################################
  
  # Test if one object is missing
  
  if(is.null(geno_off)){stop('Argument geno_off is not provided.')}
  
  if(is.null(geno_par)){stop('Argument geno_par is not provided.')}
  
  if(is.null(map)){stop('Argument map is not provided.')}
  
  if(is.null(pheno)){stop('Argument pheno is not provided.')}
  
  if(is.null(cross_ind)){stop('Argument cross_ind is not provided.')}
  
  if(is.null(par_per_cross)){stop('Argument par_per_cross is not provided.')}
  
  # test format of the marker matrices (geno_off, geno_par)
  
  if(!is.matrix(geno_off)){ stop("The geno_off argument is not a matrix.") }
  
  if(!is.character(geno_off)){
    
    stop("The geno_off argument is not a character matrix.")
    
  }
  
  if(!is.matrix(geno_par)){ stop("The geno_par argument is not a matrix.") }
  
  if(!is.character(geno_par)){
    
    stop("The geno_par argument is not a character matrix.")
    
  }
  
  # test if the marker identifiers are the same in the map and in the marker
  # matrices
  
  if(!identical(colnames(geno_off), map[, 1])){
    
    stop(paste("The marker identifier of the offspring marker matrix (geno_off)",
               "and the list of marker (map[, 1]) are not the same."))
    
  }
  
  if(!identical(colnames(geno_par), map[, 1])){
    
    stop(paste("The marker identifier of the parent marker matrix (geno_par)",
               "and the list of marker (map[, 1]) are not the same."))
    
  }
  
  # test phenotype values
  
  if(!is.matrix(pheno)){stop('The pheno argument is not a matrix.')}
  
  if(!is.numeric(pheno)){stop('The pheno argument is not a numeric matrix.')}
  
  # test if the list of genotype is the same between the offspring marker matrix
  # and the phenotypic values.
  
  if(!identical(rownames(geno_off), rownames(pheno))){
    
    stop(paste("The genotype identifiers of the offspring marker matrix",
               "(geno_off) and the list of genotypes (rownames(pheno)) are not",
               "the same."))
    
  }
  
  # length cross_ind same as list of genotypes
  
  if(length(cross_ind) != dim(geno_off)[1]){
    
    stop(paste("The cross indicator vector (cross_ind) length does not have",
               "the same length as the genotype list (dim(geno_off)[1])"))
  }
  
  # test par_per_cross
  
  if(!is.matrix(par_per_cross)){
    
    stop("The par_per_cross argument is not a matrix.")
    
  }
  
  if(!is.character(par_per_cross)){
    
    stop("The par_per_cross argument is not a character matrix.")
    
  }
  
  # remove the eventual rownames of par_per_cross
  
  if(!is.null(rownames(par_per_cross))){
    
    rownames(par_per_cross) <- NULL
    
  }
  
  if (!identical(unique(cross_ind), par_per_cross[, 1])){
    
    stop(paste("The cross identifiers used in cross_ind and",
               "in par_per_cross are different."))
    
  }
  
  # test the similarity of parents list between par_per_cross and
  # rownames(geno_par)
  
  parents <- union(par_per_cross[, 2], par_per_cross[, 3])
  
  if(sum(!(parents %in% rownames(geno_par)))>0){
    
    list.par <- paste(parents[!(parents %in% rownames(geno_par))])
    
    stop(paste("The following parents indicators:", list.par,
               "(is) are present in par_per_cross object but not in the",
               "rownames of the geno_par matrix"))
    
  }
  
  # Check the number of connected parts
  
  con_part <- design_connectivity(par_per_cross, plot_des = FALSE)
  n_con <- length(con_part)
  
  n_geno <- dim(geno_off)[1]
  n_par <- dim(geno_par)[1]
  n_cr <- dim(par_per_cross)[1]
  n_pheno <- dim(pheno)[2]
  
  # Possibility to check here the necessity to have at least 2 crosses and 3
  # different parents.
  
  ####### end check format
  
  mppData <- list(geno_off = geno_off, geno_IBS = NULL, geno_IBD = NULL,
                  geno_id = NULL, ped_mat = NULL, allele_ref = NULL,
                  geno_par = geno_par, geno_par_clu = NULL, par_clu = NULL,
                  pheno = pheno, map = map, haplo_map = NULL,
                  cross_ind = cross_ind, par_per_cross = par_per_cross,
                  type = NULL, parents = NULL, n.cr = NULL, n.par = NULL,
                  n.zigo = NULL, rem.mk = NULL, rem.gen = NULL,
                  complete = FALSE)
  
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