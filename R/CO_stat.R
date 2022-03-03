###########
# CO_stat #
###########

#' Cross-over statistics
#' 
#' Function that counts the total number of cross-over (CO) per cross, the number of cross
#' over per cross per genotype and the number of cross-over.
#'
#' @param mppData An object of class \code{mppData}
#' 
#' @return
#' 
#' \code{data.frame} with the number of cross-over per cross, average
#' number of cross-over per cross per genotype, and the average number of
#' cross-over per genotype per chromosome.
#' 
#' @author Vincent Garin
#'
#' @examples
#' 
#' # Later
#'
#' @export
#'

CO_stat <- function(mppData){
  
  cr_obj <- mppData$geno.IBD
  cr_ind <- mppData$cross.ind
  nchr = length(cr_obj$geno)
  
  cr_obj$pheno <- data.frame(Genotype = factor(paste0('G', 1:nrow(cr_obj$pheno))))
  
  N_CO <- profileGen(cross = cr_obj, bychr = TRUE, stat.type = 'xo')
  N_CO <- data.frame(cr_ind, N_CO$stat$xo)
  
  N_CO <- N_CO %>% rowwise() %>% mutate(tot_CO = sum(c_across(2:(nchr + 1)))) %>%
    group_by(cr_ind) %>% summarise(N_CO = sum(tot_CO),
                                   CO_geno = mean(tot_CO),
                                   CO_geno_chr = mean(tot_CO)/nchr)
  
  return(N_CO)
  
}