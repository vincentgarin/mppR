############
# mpp_perm #
############

#' QTL significance threshold by permutation
#'
#' Determination of an empirical null distribution of the QTL significance
#' threshold for a MPP QTL analysis using permutation test
#' (Churchill and Doerge, 1994).
#' 
#' Performs N permutations of the trait data and
#' computes each time a genome-wide QTL profile. For every run, it stores the
#' highest -log10(p-val). These values can be used to build a null distribution
#' for the QTL significance thershold. Quantile values can be determined from
#' the previous distribution. For more details about the different possible
#' models and their assumptions see \code{\link{mpp_SIM}} documentation.
#' 
#' \strong{WARNING!(1)} The computation of \code{mpp_perm()} function using mixed
#' models (all models with \code{VCOV} different than \code{"h.err"})
#' is technically possible but can be irrealistic in practice due to a reduced
#' computer power. Since a mixed model is computed at each single position it
#' can take a lot of time. From our estimation it can take between 20 to 50
#' times more time than for the linear model (HRT). We advice to compute
#' threshold with the HRT (linear) model and use it for the mixed model as well.
#' 
#' \strong{WARNING!(2)} The computation of random pedigree models
#' (\code{VCOV = "pedigree" and "ped_cr.err"}) can sometimes fail. This could be
#' due to singularities due to a strong correlation between the QTL term(s) and 
#' the polygenic term. This situation can appear in the parental model.
#' the error can also sometimes come from the \code{asreml()} function. From
#' our experience, in that case, trying to re-run the function one or two times
#' allow to obtain a result.
#'
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effects: 1) "cr" for cross-specific; 2) "par" for parental; 3) "anc"
#' for ancestral; 4) "biall" for a bi-allelic. For more details see
#' \code{\link{mpp_SIM}}. Default = "cr".
#'
#' @param VCOV \code{Character} expression defining the type of variance
#' covariance structure used: 1) "h.err" for an homogeneous variance residual term
#' (HRT) linear model; 2) "h.err.as" for a HRT model fitted by REML using
#' \code{ASReml-R}; 3) "cr.err" for a cross-specific variance residual terms
#' (CSRT) model; 4) "pedigree" for a random pedigree term and HRT model;
#' and 5) "ped_cr.err" for random pedigree and CSRT model.
#' For more details see \code{\link{mpp_SIM}}. Default = "h.err".
#' 
#' @param N Number of permutations. Default = 1000.
#' 
#' @param q.val Single \code{numeric} valu or vector of desired quantiles from
#' the null distribution. Default = 0.95.
#' 
#' @param verbose \code{Logical} value indicating if progression of the function
#' should be printed. It will not affect the printing of the other functions
#' called by \code{mpp_perm()}, especially the printing of \code{asreml()}.
#' Default = TRUE.
#' 
#' @param n.cores \code{Numeric}. Specify here the number of cores you like to
#' use. Default = 1.
#'
#' @return Return:
#'
#' \code{List} with the following object:
#'
#' \item{max.pval }{Vector of the highest genome-wide -log10(p-values).}
#'
#' \item{q.val }{Quantile values from the QTL significance threshold null
#' distribution.}
#'
#' \item{seed}{\code{Numeric} vector of random generated seed values for each
#' permutation.}
#'
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mpp_SIM}}
#'
#' @references
#'
#' Churchill, G. A., & Doerge, R. W. (1994). Empirical threshold values for
#' quantitative trait mapping. Genetics, 138(3), 963-971.
#'
#' @examples
#' 
#' \dontrun{
#'
#' data(mppData)
#' 
#' Perm <- mpp_perm(mppData = mppData, Q.eff = "cr", VCOV = "h.err",
#'                        N = 100)
#' 
#' }
#'
#'
#' @export
#'


mpp_perm <- function(mppData, trait = 1, Q.eff = 'cr', VCOV = "h.err", N = 1000,
                     q.val = 0.95, verbose = TRUE, n.cores = 1) {
  
  # 1. Check data format and arguments
  ####################################
  
  check.model.comp(mppData = mppData, trait = trait, Q.eff = Q.eff, VCOV = VCOV,
                   n.cores = n.cores, fct = "perm")
  
  
  # 2. Form required elements for the analysis
  ############################################
  
  ### 2.1 trait values
  
  t_val <- sel_trait(mppData = mppData, trait = trait)
  
  ### 2.2 inverse of the pedigree matrix
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.3 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.4 create space to store the results
  
  max.pval <- numeric(N)
  seed <- numeric(N)
  vect.pos <- 1:dim(mppData$map)[1]
  
  if (N >= 100) { reference <- reference.count(N = N, l = 10) ;count <- 1 }
  
  ### 2.5 Optional cluster
  
  if(n.cores > 1){
    
    parallel <- TRUE
    cluster <- makeCluster(n.cores)
    
  } else {
    
    parallel <- FALSE
    cluster <- NULL
    
  }
  
  # 3. Run the permutations
  #########################
  
  for (i in 1:N) {
    
    seed[i] <- round(runif(n = 1, min = 1, max = 1e+06))
    set.seed(seed[i])
    
    ### 3.1 trait permutation (within cross)
    
    perm_cross <- function(x) x[sample(1:length(x))]
    
    cross.ind.fac <- factor(x = mppData$cross.ind,
                            levels = unique(mppData$cross.ind))
    
    t_val_i <- unlist(tapply(t_val, cross.ind.fac, perm_cross))
    
    ### 3.2 genome scan 
    
    if (parallel) {
      
      perm.i <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelPerm, 
                          mppData = mppData, trait = t_val_i,
                          cross.mat = cross.mat, Q.eff = Q.eff, VCOV = VCOV)
      
    } else {
      
      perm.i <- lapply(X = vect.pos, FUN = QTLModelPerm,  mppData = mppData,
                       trait = t_val_i, cross.mat = cross.mat, Q.eff = Q.eff,
                       VCOV = VCOV)
      
    }
    
    max.pval[i] <- max(unlist(perm.i), na.rm = TRUE)
    
    # flag the progresses of the process
    
    if(verbose){
      
      if (N >= 100) {
        
        if (i >= reference[count, 1]) {
          
          cat(reference[count, 2], "%")
          cat("\n")
          count <- count + 1
          
        }
        
      }
      
    }
    
    
  } # end permutations
  
  if(n.cores > 1){stopCluster(cluster)}
  
  # 4. return.results
  ###################
  
  # quantile values
  
  q.val <- quantile(max.pval, q.val)
  
  if(verbose){
    
    hist(max.pval, nclass = 20, main = "histogram of the maxium -log10(pvalue)")
    
    cat("\n")
    cat("Desired quantiles:\n")
    cat("\n")
    
    print(q.val)
    
  }
  
  results <- list(max.pval = max.pval, q.val = q.val, seed = seed)
  
  return(results)
  
}