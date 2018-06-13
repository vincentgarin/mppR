#################
# formPedMatInv #
#################

# function that form the inverse of the pedigree matrix for asreml function and
# export it in the global environment.

formPedMatInv <- function(mppData, VCOV){
  
  # if((VCOV == "pedigree") || (VCOV == "ped_cr.err")){
  #   
  #   assign("ped.mat.inv", asreml::asreml.Ainverse(mppData$ped.mat[, 2:4]),
  #          envir = .GlobalEnv)
  #   
  #   # extract the sparse inverse
  #   
  #   assign("ped.mat.inv", ped.mat.inv$ginv, envir = .GlobalEnv)
  #   
  # }
  
}