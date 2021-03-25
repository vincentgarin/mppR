###################
# replace_prob_mk #
###################

replace_prob_mk <- function(geno.off, geno.par, prob.mk, reference_allele){
  
  for(i in 1:length(prob.mk)){}
  
  ref_i <- reference_allele[reference_allele[, 1] == prob.mk[i], 2]
  ref_i <- unlist(strsplit(x = gsub('/', '', ref_i), split = ''))
  mk_sc_i <- c(paste0(ref_i[1], ref_i[1]),
               paste0(ref_i[2], ref_i[2]),
               paste0(ref_i[1], ref_i[2]),
               paste0(ref_i[2], ref_i[1]))
  
  geno.off[!(geno.off[, prob.mk[i]] %in% mk_sc_i), prob.mk[i]] <- NA
  
  geno.par[!(geno.par[, prob.mk[i]] %in% mk_sc_i), prob.mk[i]] <- NA
  
  return(list(geno.off = geno.off, geno.par = geno.par))
  
}