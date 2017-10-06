################
# QC_matchGeno #
################

#' Match genotypes in the marker matrix and the phenotype data
#' 
#' Determine which genotypes are common between the genotype marker matrix
#' and the phenotypic data. Return the marker matrix and the phenotypic values
#' that correspond to these genotypes.
#' 
#' 
#' @param mk.mat Marker score \code{matrix} or \code{data.frame} with genotypes
#' as row and markers as column. \strong{Genotype indicators must be the row
#' names of the matrix.}
#' 
#' @param pheno Matrix or data frame representing phenotypic data
#' \strong{with genotype indicators as row names}.
#' 
#' @return Return:
#' 
#' \item{new.mk.mat}{Marker \code{matrix} with only the common genotypes.}
#' 
#' \item{new.pheno}{Phenotypic data of the common genotypes.}
#' 
#' \item{mk.mat.only}{Markers present in the marker matrix but not in
#' the phenotypic data.}
#' 
#' \item{pheno.only}{Markers present in the phenotypic data but not in
#' the marker matrix.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{QC_matchMarker}}
#' 
#' @examples
#' 
#' data(USNAM_pheno)
#' data(USNAM_geno)
#' 
#' rownames(USNAM_pheno) <- USNAM_pheno$genotypes
#' 
#' # The 6 parents of the crosses were not phenotyped
#' match <- QC_matchGeno(mk.mat = USNAM_geno, pheno = USNAM_pheno)
#' 
#' # Get the consensus phenotype data and genotype matrix
#' new.mk.mat <- match$new.mk.mat
#' new.pheno <- match$new.pheno
#' rm(match)
#' 
#' @export
#' 


QC_matchGeno <- function(mk.mat, pheno){
  
  # identify common marker in marker matrix and in map (intersection A n B)
  
  inter.geno.pheno <- intersect(rownames(mk.mat), rownames(pheno))
  
  # genotypes only in the genotype matrix (A \B)
  
  geno.only <- setdiff(rownames(mk.mat), rownames(pheno))
  
  # genotypes only in the phenotype file (B \A)
  
  pheno.only <- setdiff(rownames(pheno), rownames(mk.mat))
  
  # subset the marker matrix and the phenotype list
  
  new.mk.mat <- subset(x = mk.mat, subset = rownames(mk.mat) %in% inter.geno.pheno,
                       drop = FALSE)
  
  new.pheno <- subset(x = pheno, subset = rownames(pheno) %in% inter.geno.pheno,
                      drop = FALSE)
  
  # make sure the order is the same in genotype and phenotype
  
  new.pheno <- new.pheno[rownames(new.mk.mat), ]
  
  return(list(new.mk.mat=new.mk.mat, new.pheno=new.pheno, mk.mat.only=geno.only,
              pheno.only=pheno.only))
  
}