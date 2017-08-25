############
# mpp_perm #
############

#' QTL significance threshold by permutation
#'
#' Determination of the null distribution of the QTL significance threshold for
#' a MPP QTL analysis using permutation test (Churchill and Doerge, 1994).
#' For more details about the different possible models and their assumptions
#' see \code{\link{mpp_SIM}} documentation.
#' 
#' Performs N permutations of the trait data and
#' computes each time a genome-wide QTL profile. For every run, it stores the
#' highest -log10(p-val). These values can be used to build a null distribution
#' for the QTL significance thershold. Quantile values can be determined from
#' the previous distribution.
#' 
#' \strong{WARNING!(1)} The computation of \code{mpp_perm()} function using mixed
#' models (all models with \code{VCOV} different than \code{"h.err"})
#' is technically possible but can be irrealistic
#' in practice due to a reduced computer power. Since a mixed model is computed at
#' each single position it can take a lot of time. From our estimation it can take
#' between 20 to 50 times more time than for linear models. We advice to compute
#' threshold with the HRT (linear) model and use it for the mixed model as well.
#' 
#' \strong{WARNING!(2)} The estimation of the random pedigree models
#' (\code{VCOV = "pedigree" and "ped_cr.err"}) can be unstable. Sometimes the
#' \code{asreml()} function fails to produce a results and returns the following
#' message: \strong{\code{GIV matrix not positive definite: Singular pivots}}.
#' So far we were not able to identify the reason of this problem and to
#' reproduce this error because it seems to happen randomly. The consequence of
#' this is that part of the results of the CV procedure can not be produced.
#' These will be replaced by 0 values.
#' 
#' 
#'
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. For more details see \code{\link{mpp_SIM}}. Default = "cr".
#'
#' @param par.clu Required argument for the ancesral model \code{(Q.eff = "anc")}.
#' \code{interger matrix} representing the results of a parents genotypes
#' clustering. The columns represent the parental lines and the rows
#' the different markers or in between positions. \strong{The columns names must
#' be the same as the parents list of the mppData object. The rownames must be
#' the same as the map marker list of the mppData object.} At a particular
#' position, parents with the same value are assumed to inherit from the same
#' ancestor. for more details, see \code{\link{USNAM_parClu}} and
#' \code{\link{parent_cluster}}. Default = NULL.
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
#' @param parallel \code{Logical} value specifying if the function should be
#' executed in parallel on multiple cores. To run function in parallel user must
#' provide cluster in the \code{cluster} argument. \strong{Parallelization is
#' only available for HRT (linear) models \code{VCOV = "h.err"}}.
#' Default = FALSE.
#'
#' @param cluster Cluster object obtained with the function \code{makeCluster()}
#' from the parallel package. Default = NULL.
#' 
#' @param verbose \code{Logical} value indicating if progression of the function
#' should be printed. It will not affect the printing of the other functions
#' called by \code{mpp_perm()}, especially the printing of \code{asreml()}.
#' Default = TRUE.
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
#' @seealso \code{\link{mppData_form}}, \code{\link{mpp_SIM}},
#' \code{\link{parent_cluster}},
#' \code{\link{USNAM_parClu}} 
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
#' data(USNAM_mppData)
#' 
#' PermTestCr <- mpp_perm(mppData = USNAM_mppData, Q.eff = "cr", VCOV = "h.err",
#'                        N = 100)
#' 
#' # Using multiple core
#' 
#' library(parallel)
#' n.cores <- detectCores()
#' cluster <- makeCluster((n.cores-1))
#' 
#' PermTestCr <- mpp_perm(mppData = USNAM_mppData, Q.eff = "cr", VCOV = "h.err",
#'                        N = 100, parallel = TRUE, cluster = cluster)
#' 
#' stopCluster(cl = cluster)
#' 
#' }
#'
#'
#' @export
#'


mpp_perm <- function(mppData, Q.eff, par.clu = NULL, VCOV = "h.err", N = 1000, 
                     q.val = 0.95, parallel = FALSE, cluster = NULL,
                     verbose = TRUE) {
  
  # 1. Check data format and arguments
  ####################################
  
  check.model.comp(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
                   par.clu = par.clu, parallel = parallel, cluster = cluster,
                   fct = "perm")
  
  
  # 2. Form required elements for the analysis
  ############################################
  
  ### 2.1 inverse of the pedigree matrix
  
  formPedMatInv(mppData = mppData, VCOV = VCOV)
  
  ### 2.2 cross matrix (cross intercept)
  
  cross.mat <- IncMat_cross(cross.ind = mppData$cross.ind)
  
  ### 2.3 parent matrix
  
  parent.mat <- IncMat_parent(mppData)
  
  ### 2.4 modify the par.clu object order parents columns and replace monomorphic
  
  if (Q.eff == "anc") {
    
    check <- parent_clusterCheck(par.clu = par.clu)
    par.clu <- check$par.clu[, mppData$parents] # order parents columns
    
  } else {par.clu <- NULL}
  
  
  ### 2.5 create space to store the results
  
  max.pval <- numeric(N)
  seed <- numeric(N)
  vect.pos <- 1:dim(mppData$map)[1]
  
  if (N >= 100) { reference <- reference.count(N = N, l = 10) ;count <- 1 }
  
  
  # 3. Run the permutations
  #########################
  
  for (i in 1:N) {
    
    seed[i] <- round(runif(n = 1, min = 1, max = 1e+06))
    set.seed(seed[i])
    
    ### 3.1 trait permutation (within cross)
    
    perm_cross <- function(x) x[sample(1:length(x))]
    
    cross.ind.fac <- factor(x = mppData$cross.ind,
                            levels = unique(mppData$cross.ind))
    
    mppData$trait[, 1] <- unlist(tapply(mppData$trait[, 1], cross.ind.fac,
                                         perm_cross))
    
    ### 3.2 genome scan 
    
    if (parallel) {
      
      perm.i <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelPerm, 
                       mppData = mppData, cross.mat = cross.mat,
                       par.mat = parent.mat, Q.eff = Q.eff,
                       par.clu = par.clu, VCOV = VCOV)
      
    } else {
      
      perm.i <- lapply(X = vect.pos, FUN = QTLModelPerm, 
                      mppData = mppData, cross.mat = cross.mat,
                      par.mat = parent.mat, Q.eff = Q.eff,
                      par.clu = par.clu, VCOV = VCOV)
        
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