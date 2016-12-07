###########
# mpp_SIM #
###########

#' MPP Simple Interval Maping
#' 
#' Computes single QTL models with different possible assumptions concerning
#' the number of alleles at the QTL position and the variance covariance
#' structure (VCOV) of the model. The function returns a -log10(p-value) QTL
#' profile. 
#' 
#' The implemented models vary according to the number of alleles assumed at the
#' QTL position (and their origin) and their variance covariance structure.
#' In both case four assumptions are possible giving a grid of 16 different
#' models.
#' 
#' Concerning the type of QTL effect, the first option a cross-specific QTL
#' effects model (\code{Q.eff = "cr"}). In this model, the QTL effects are
#' assumed to be nested within cross which leads to the estimation of one
#' parameter per cross. Computed with a homogoenous residual variance term
#' (\code{VCOV = "h.err"}), the cross-specific model corresponds to the
#' disconnected described in Blanc et al. 2006.
#' 
#' A second possibility is the parental model (\code{Q.eff = "par"}). The
#' parental model assumes 1 QTL effect (allele) per parent that are independent
#' from the genetic background. This means that QTL coming form parent i has the
#' same effect in all crosses where this parent is used. This model is suppose to
#' allow to have better estimates of the QTL due to larger sample size due to
#' shared parents. If the number of parents is lower than the number of cross,
#' it is also expected to be more powerfull due to a reduced number of QTL
#' parameters to estimate assuming that the hypothesis of similar parental effects
#' through crosses holds. Calculated with HRT assumption, the parental model
#' correspond to the connected model presented in Blanc et al. (2006).
#' 
#' The third type of model is the ancestral model (code{Q.eff = "anc"}). This
#' model tries to use genetic relatedness that could exist between parents.
#' Indeed, the parental model assumes that parent are independent which is not
#' the case. Using genetic relatedness between the parents, it is possible group
#' these parents into a reduced number of ancestral cluster. Parents belonging
#' to the same ancestral group are assumed to transmit the same allele
#' (Jansen et al. 2003; Leroux et al. 2014). Clustering of parental line can
#' be performed using the \code{\link{parent_cluster}} function which calls
#' functions of the R package
#' \code{clustahplo}. This information allows to modify the IBD relationship
#' between the parents. The ancestral model estimate therefore one QTL effect
#' per ancestral class. Once again, the theoretical expectation is a gain of
#' QTL detection power by the reduction of the number of parameters to estimate.
#' The HRT ancestral model correspond to the linkage desequilibrium
#' linkage analysis (LDLA) models used by Bardol et al. (2013) or
#' Giraud et al. (2014).
#' 
#' The final posibility is the bi-allelic model (\code{Q.eff = "biall"}).
#' Bi-allelic genetic predictor are a single vector with value 0, 1 or 2
#' corresponding to the number of allele copy of the most frequent SNP allele.
#' Relatedness between line is therefore defined via identical by state (IBS)
#' measurement. This model corresponds to models used for association mapping.
#' For example, if (\code{VCOV = "h.err"}), it is similar to model B in
#' Wurschum et al. (2012) or association mapping model in Liu et al. (2012).
#' 
#' The second dimension relates to the form of the variance covariance structure.
#' The first assumption is an homogeneous variance residual term model (HRT). It
#' corresponds to the classical linear model assumption of independence of the
#' residual terms. It can be calculated using the \code{lm()} function
#' (\code{VCOV = "h.err"}) or by REML using \code{asreml()}
#' (\code{VCOV = "h.err.as"}).
#' 
#' The second assumption is a cross-specific variance residual terms model
#' (CSRT) where one residual variance component is estimated per cross 
#' (\code{VCOV = "cr.err"}). This VCOV allows to take into consideration
#' potential heterogeneity of residual variance that could exist between crosses.
#' Such a difference can be due to different levels of polygenic effect
#' (undetected QTLs).
#' 
#' The third and fourth VCOV model the polygenic effect directly by adding a
#' random pedigree term to the model. The variance covariance structure
#' associated with this random term is the matrix of inbreeding coefficients
#' based on pedigree information (A). It is calculated using the function
#' \code{asreml.Ainverse()}. This function uses the method developped by
#' Meuwissen and Luo (1992). The third possibility is to use the random pedigree
#' model with HRT (\code{VCOV = "pedigree"}). The fourth option includes the
#' random pedigree term plus CSRT (\code{VCOV = "ped_cr.err"}). All model except
#' (\code{VCOV = "h.err"}) are fitted using REML via the \code{asreml()} function
#' from the \code{ASReml-R} package (Butler et al., 2009).
#' 
#' \strong{WARNING!} The estimation of the random pedigree models
#' (\code{VCOV = "pedigree" and "ped_cr.err"}) can be unstable. Sometimes the
#' \code{asreml()} function fails to produce a results and returns the following
#' message: \strong{\code{GIV matrix not positive definite: Singular pivots}}.
#' So far we were not able to identify the reason of this problem and to
#' reproduce this error because it seems to happen randomly. From our
#' experience, trying to re-run the function one or two times should allow
#' to obtain a result.
#' 
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#'
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects; 4) "biall" for a bi-allelic
#' effects. Default = "cr".
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
#' Default = "h.err".
#' 
#' @param est.gen.eff \code{Logical} value. if \code{est.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot_genEffects}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' \strong{This functionality is ony available for the cross-specific,
#' parental and ancestral models.}
#' Default value = FALSE.
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
#'   
#' @return Return:
#' 
#' \item{SIM }{\code{Data.frame} of class \code{QTLprof}. with five columns :
#' 1) QTL marker or in between position names; 2) chromosomes;
#' 3) Interger position indicators on the chromosome;
#' 4) positions in centi-Morgan; and 5) -log10(p-val). And if
#' \code{est.gen.eff = TRUE}, p-values of the cross or parental QTL effects.}
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{mppData_form}}, \code{\link{parent_cluster}},
#' \code{\link{plot_genEffects}},
#' \code{\link{USNAM_parClu}}
#' 
#' @references
#' 
#' Butler D. G., Cullis B. R., Gilmour A. R., Gogel B. J. (2009) ASReml-R
#' reference manual.
#' 
#' Bardol, N., Ventelon, M., Mangin, B., Jasson, S., Loywick, V., Couton, F., ...
#' & Moreau, L. (2013). Combined linkage and linkage disequilibrium QTL mapping 
#' in multiple families of maize (Zea mays L.) line crosses highlights
#' complementarities between models based on parental haplotype and single locus
#' polymorphism. Theoretical and applied genetics, 126(11), 2717-2736.
#' 
#' Blanc, G., Charcosset, A., Mangin, B., Gallais, A., & Moreau, L. (2006).
#' Connected populations for detecting quantitative trait loci and testing for
#' epistasis: an application in maize. Theoretical and Applied Genetics,
#' 113(2), 206-224. 
#' 
#' Giraud, H., Lehermeier, C., Bauer, E., Falque, M., Segura, V., Bauland,
#' C., ... & Moreau, L. (2014). Linkage Disequilibrium with Linkage Analysis
#' of Multiline Crosses Reveals Different Multiallelic QTL for Hybrid
#' Performance in the Flint and Dent Heterotic Groups of Maize. Genetics,
#' 198(4), 1717-1734.
#' 
#' Jansen, R. C., Jannink, J. L., & Beavis, W. D. (2003). Mapping quantitative
#' trait loci in plant breeding populations. Crop Science, 43(3), 829-834.
#' 
#' Leroux, D., Rahmani, A., Jasson, S., Ventelon, M., Louis, F., Moreau, L.,
#' & Mangin, B. (2014). Clusthaplo: a plug-in for MCQTL to enhance QTL detection
#' using ancestral alleles in multi-cross design. Theoretical and Applied
#' Genetics, 127(4), 921-933.
#' 
#' Liu, W., Reif, J. C., Ranc, N., Della Porta, G., & Wurschum, T. (2012).
#' Comparison of biometrical approaches for QTL detection in multiple
#' segregating families. Theoretical and Applied Genetics, 125(5), 987-998.
#' 
#' Meuwissen T and Luo, Z. (1992). Computing inbreeding coefficients in large
#' populations. Genetics Selection Evolution, 24(4), 305-313.
#' 
#' Wurschum, T., Liu, W., Gowda, M., Maurer, H. P., Fischer, S., Schechert, A.,
#' & Reif, J. C. (2012). Comparison of biometrical models for joint linkage
#' association mapping. Heredity, 108(3), 332-340. 
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' 
#' # Cross-specific model
#' ######################
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData, Q.eff = "cr", VCOV = "h.err",
#' est.gen.eff = TRUE)
#' 
#' plot_QTLprof(Qprof = SIM)  
#' plot_genEffects(Qprof = SIM)
#' 
#' \dontrun{
#' 
#' # Using multiple cores
#' 
#' library(parallel)
#' n.cores <- detectCores()
#' cluster <- makeCluster((n.cores-1))
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData, Q.eff = "cr", VCOV = "h.err",
#'                parallel = TRUE, cluster = cluster)
#' 
#' stopCluster(cl = cluster)
#' 
#' }
#' 
#' # Bi-allelic model
#' ##################
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData_bi, Q.eff = "biall", VCOV = "h.err")
#' 
#' plot_QTLprof(Qprof = SIM, type = "h")
#' 
#' @export
#' 


mpp_SIM <- function(mppData, Q.eff = "cr", par.clu = NULL, VCOV = "h.err", 
                    est.gen.eff = FALSE, parallel = FALSE, cluster = NULL) {
  
  # 1. Check data format and arguments
  ####################################
  
  check.model.comp(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
                   par.clu = par.clu, est.gen.eff = est.gen.eff,
                   parallel = parallel, cluster = cluster,
                   fct = "SIM")
  
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
  
  vect.pos <- 1:dim(mppData$map)[1]
  
  # 3. computation of the SIM profile (genome scan)
  #################################################
  
  if (parallel) {
    
    log.pval <- parLapply(cl = cluster, X = vect.pos, fun = QTLModelSIM,
                        mppData = mppData, cross.mat = cross.mat,
                        par.mat = parent.mat, Q.eff = Q.eff,
                        par.clu = par.clu, VCOV = VCOV,
                        est.gen.eff = est.gen.eff)
    
  } else {
    
    log.pval <- lapply(X = vect.pos, FUN = QTLModelSIM,
                     mppData = mppData, cross.mat = cross.mat,
                     par.mat = parent.mat, Q.eff = Q.eff,
                     par.clu = par.clu, VCOV = VCOV,
                     est.gen.eff = est.gen.eff)
    
  }
  
  log.pval <- t(data.frame(log.pval))
  log.pval[, 1] <- check.inf(x = log.pval[, 1]) # check if there are -/+ Inf value
  log.pval[is.na(log.pval[, 1]), 1] <- 0
  
  # 4. form the results
  #####################
  
  SIM <- data.frame(mppData$map, log.pval)
  
  
  if(est.gen.eff){
    
    if(Q.eff=="cr"){ Qeff_names <- unique(mppData$cross.ind)
      
      } else { Qeff_names <- mppData$parents }
    
    colnames(SIM)[5:dim(SIM)[2]] <- c("log10pval", Qeff_names)
    
  } else {colnames(SIM)[5] <- "log10pval"} 
  
  class(SIM) <- c("QTLprof", "data.frame")
  
  ### 4.1: Verify the positions for which model could not be computed
  
  if(sum(SIM$log10pval == 0) > 0) {
    
    if (sum(SIM$log10pval) == 0){
      
      message("The computation of the QTL models failled for all positions.
              This could be due to problem in asreml function.")
      
    } else {
      
      list.pos <- mppData$map[(SIM$log10pval == 0), 1]
      
      end.mess <- ". This could be due to singularities or function issue."
      
      text <- paste("The computation of the QTL model failed for the following positions: ",
                    paste(list.pos, collapse = ", "), end.mess)
      
      message(text)
      
    }
    
  }
  
  return(SIM)
  
}