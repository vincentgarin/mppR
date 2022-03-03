#################
# geno_ABH_prob #
#################

#' Probability of the genotype
#' 
#' Function that calculate the probability of parental origin (AA), (AB) or (BB)
#' genotypes.
#'
#' @param mppData An object of class \code{mppData}
#' 
#' @return
#' 
#' Vector of parental orignin genotype probability.
#' 
#' @author Vincent Garin
#'
#' @examples
#' 
#' # Later
#'
#' @export
#'

geno_ABH_prob <- function(mppData){
  
  cr_obj <- mppData$geno.IBD
  
  tot_1 <- 0
  tot_2 <- 0
  tot_3 <- 0
  
  for(i in 1:length(cr_obj$geno)){
    
    d_i <- cr_obj$geno[[i]]$data
    f_i <- table(c(d_i))
    
    if(!is.na(f_i['1'])){tot_1 <- tot_1 + f_i['1']}
    if(!is.na(f_i['1'])){tot_2 <- tot_2 + f_i['2']}
    if(!is.na(f_i['1'])){tot_3 <- tot_3 + f_i['3']}
    
  }
  
  f_ABH <- c(tot_1, tot_2, tot_3)/sum(c(tot_1, tot_2, tot_3))
  names(f_ABH) <- c('P(AA)', 'P(AB)', 'P(BB)')
  
  return(f_ABH)
  
}