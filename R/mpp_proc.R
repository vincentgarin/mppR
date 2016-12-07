############
# mpp_proc #
############

#' MPP QTL analysis
#' 
#' Multi-parent QTL analysis based on models with different possible assumptions
#' concerning the number of alleles at the QTL position and the variance
#' covariance structure (VCOV) of the model. For more details about the different
#' models, see documentation of the function \code{\link{mpp_SIM}}.
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
#' The procedure is the following:
#' 
#' \enumerate{
#' 
#' \item{Simple interval mapping (SIM) to select cofactor
#' (\code{\link{mpp_SIM}}).}
#' 
#' \item{Composite interval mapping (CIM) with selected cofactors
#' (\code{\link{mpp_CIM}}).}
#' 
#' \item{Optional backward elimination on the list of QTL
#' candidates (\code{backward = TRUE}) (\code{\link{mpp_BackElim}}).}
#' 
#' \item{Computation of the QTL genetic effects (\code{\link{QTL_genEffects}})
#' and proportion of the phenotypic variation explained by the QTLs (R squared)
#' (\code{\link{QTL_R2}}).}
#' 
#' \item{Optional QTL confidence interval computation from a CIM- profile
#' (excluding cofactors on the scanned chromosome) (\code{argument CI=TRUE})
#' (\code{\link{QTL_CI}}).}
#' 
#' }
#' 
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP".
#' 
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#' 
#' @param mppData An object of class \code{mppData}. See
#' \code{\link{mppData_form}} for details.
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
#' @param est.gen.eff \code{Logical} value. if \code{est.gen.eff = TRUE},
#' the function will save the decomposed genetic effects per cross/parent.
#' These results can be ploted with the function \code{\link{plot_genEffects}}
#' to visualize a genome-wide decomposition of the genetic effects.
#' \strong{This functionality is ony available for the cross-specific,
#' parental and ancestral models.}
#' Default value = FALSE.
#' 
#' @param thre.cof \code{Numeric} value representing the -log10(p-value)
#' threshold
#' above which a position can be peaked as a cofactor. Significance threshold
#' values can be obtained by permutation using \code{\link{mpp_perm}} function.
#' Default = 3.
#' 
#' @param win.cof \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected cofactors. Default = 20.
#' 
#' @param N.cim \code{Numeric} value specifying the number of time the CIM
#' analysis is repeated. Default = 1.
#' 
#' @param window \code{Numeric} distance on the left an right of a cofactor
#' position where it is not included in the model. Default = 20.
#' 
#' @param thre.QTL \code{Numeric} value representing the -log10(p-value)
#' threshold above which a position can be selected as QTL. Significance
#' threshold values can be obtained by permutation using \code{\link{mpp_perm}}
#' function. Default = 3.
#' 
#' @param win.QTL \code{Numeric} value in centi-Morgan representing the minimum
#' distance between two selected QTL. Default = 20.
#' 
#' @param backward \code{Logical} value. If \code{backward = TRUE},
#' the function performs a backward elimination on the list of selected QTLs.
#' Default = TRUE.
#' 
#' @param alpha.bk \code{Numeric} value indicating the significance level for
#' the backward elimination. Terms with p-values above this value will
#' iteratively be removed. Default = 0.05.
#' 
#' @param LR.R2 \code{Logical} value. If \code{LR.R2 = TRUE}, the likelihood
#' ratio R squared will be used if \code{VCOV = "cr.err"}. In all other
#' situations R squared from a linear model will be computed. For more details
#' see \code{\link{QTL_R2}}. Default = TRUE.
#' 
#' @param const \code{Character} expression indicating the type of contstraint
#' used for the estimation of the genetic effect of the HRT (linear) model
#' (\code{VCOV = "h.err"}) for the parental or ancestral models
#' (\code{Q.eff = "par" or "anc"}). If const = "set.0", the value of the
#' parent or the ancestral class containing the parent specified in argument
#' (\code{par.ref}) is set to 0. If const = "sum.0", the sum of
#' the parental (ancestral) effects is constrained to sum to 0. For more details
#' about estimation of genetic effects see \code{\link{QTL_genEffects}}.
#' Default = "set.0".
#' 
#' @param par.ref \code{Character} expression indicating the parent or the
#' ancestral class containing the specified parent that will be set as reference
#' for the computation of the genetic effect. By default the function will use
#' the first parent of the \code{mppData$parents} list. Default = NULL.
#' 
#' @param CI \code{Logical} value. If \code{CI = TRUE}, the function will
#' compute a -log10(pval) drop confidence interval for each QTL after
#' calculating a CIM- profile (without cofactors on the scanned chromosome).
#' Default = FALSE.
#' 
#' @param drop \code{Numeric} -log10(p-value) drop value at the limits of the
#' interval. Default = 1.5.
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
#' @param silence.print \code{Logical} value specifying if the printing of the
#' \code{mpp_proc()} function must be silenced. It will not
#' affect the printing of the other functions called by \code{mpp_proc()},
#' especially the printing of \code{asreml()}. Default = FALSE.
#'
#' @param output.loc Path where a folder will be created to save the results.
#' By default the function uses the current working directory.
#' 
#' 
#' @return Return:
#' 
#' List containing the following items:
#' 
#' \item{n.QTL}{Number of detected QTLs}
#' 
#' \item{cofactors}{\code{Data.frame} with cofactors positions.}
#' 
#' \item{QTL}{\code{Data.frame} with QTL positions.}
#' 
#' \item{R2}{\code{List} containing R squared statistics of the QTL effects.
#' for details see \code{\link{QTL_R2}} output section.}
#' 
#' \item{QTL.effects}{\code{List} of QTLs genetic effects. For details see
#' \code{\link{QTL_genEffects}} output section.}
#' 
#' \item{QTL.CI}{If \code{CI = TRUE}, confidence interval information of
#' the QTLs.}
#' 
#' Some output files are also saved at the specified location
#' (\code{output.loc}):
#' 
#' \enumerate{
#' 
#' \item{A QTL report (QTL_REPORT.txt) with: 1) the number of detected QTLs;
#' 2) the global R squared statistics; 3) for each QTL, position information
#' (plus confidence interval if \code{CI = TRUE}) and estimated QTL genetic
#' effects per cross or parents (for details see \code{\link{QTL_genEffects}}).}
#' 
#' \item{The SIM and CIM results in a text file (SIM.txt, CIM.txt).}
#' 
#' \item{The list of cofactors (cofactors.txt).}
#' 
#' \item{The list of QTL (QTL.txt).}
#' 
#' \item{The QTL R squared statistics (QTL_R2.txt) (for details see
#' \code{\link{QTL_R2}}).}
#' 
#' \item{If \code{CI = TRUE}, the QTL confidence intervals (QTL_CI.txt).}
#' 
#' \item{General results of the QTL detection process: number of QTLs and
#' global adjusted and non-adjusted R squared statistics (QTL_genResults.txt).}
#' 
#' \item{The plot of the CIM profile (QTL_profile.pdf) with dotted vertical
#' lines representing the cofactors positions. If \code{est.gen.eff = TRUE},
#' plot of the genetic effects per cross or parents (gen_eff.pdf) with dashed
#' lines representing the QTL positions. For more details see
#' \code{\link{plot_QTLprof}} and \code{\link{plot_genEffects}}.}
#' 
#' }
#'                   
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{mppData_form}},
#' \code{\link{mpp_BackElim}},
#' \code{\link{mpp_CIM}},
#' \code{\link{mpp_perm}},
#' \code{\link{mpp_SIM}},
#' \code{\link{parent_cluster}},
#' \code{\link{plot_genEffects}},
#' \code{\link{plot_QTLprof}},
#' \code{\link{QTL_CI}},
#' \code{\link{QTL_genEffects}},
#' \code{\link{QTL_R2}},
#' \code{\link{USNAM_parClu}}
#' 
#' @examples
#'  
#' \dontrun{
#' 
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' 
#' # Detect the number of cores and build cluster
#' 
#' library(parallel)
#' n.cores <- detectCores()
#' cluster <- makeCluster(n.cores-1)
#' 
#' # Specify a location where your results will be saved
#' my.loc <- "C:/.../..."
#' 
#' # Cross-specific model
#' 
#' USNAM_cr <- mpp_proc(pop.name = "USNAM", trait.name = "ULA",
#'                      mppData = USNAM_mppData, est.gen.eff = TRUE, CI = TRUE,
#'                      parallel = TRUE, cluster = cluster, output.loc = my.loc)
#' 
#' # Bi-allelic model
#' 
#' USNAM_biall <- mpp_proc(pop.name = "USNAM", trait.name = "ULA",
#'                      mppData = USNAM_mppData_bi,Q.eff = "biall",
#'                      output.loc = my.loc)
#' 
#' }
#' 
#' 
#' @export
#'


mpp_proc <- function(pop.name = "MPP", trait.name = "trait1", mppData,
                     Q.eff = "cr", par.clu = NULL, VCOV = "h.err",
                     est.gen.eff = FALSE, thre.cof = 3, win.cof = 20,
                     N.cim = 1, window = 20, thre.QTL = 3, win.QTL = 20,
                     backward = TRUE, alpha.bk = 0.05, LR.R2 = TRUE, 
                     const = "set.0", par.ref = NULL, CI = FALSE,
                     drop = 1.5, parallel = FALSE, cluster = NULL,
                     silence.print = FALSE, output.loc = getwd()) {
  
  
# 1. Check the validity of the parameters that have been introduced
###################################################################

if(is.null(par.ref)){ par.ref <- mppData$parents[1]}

check.mpp.proc(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
               par.clu = par.clu, est.gen.eff = est.gen.eff,
               parallel = parallel, cluster = cluster, par.ref = par.ref,
               const = const, output.loc = output.loc)


# 2. Create a directory to store the results
############################################

# create a directory to store the results of the QTL analysis

end.char <- substr(output.loc, nchar(output.loc), nchar(output.loc))

if(end.char == "/"){

  folder.loc <- paste0(output.loc, paste("QTLan", pop.name, trait.name, Q.eff,
                                         VCOV, sep = "_"))

} else {

  folder.loc <- paste0(output.loc, "/", paste("QTLan", pop.name, trait.name,
                                              Q.eff, VCOV, sep = "_"))

}

dir.create(folder.loc)


# 3. Cofactors selection - SIM
##############################

if(!silence.print){

  cat("\n")
  cat("Cofactors selection - SIM")
  cat("\n")
  cat("\n")

}

SIM <- mpp_SIM(mppData = mppData, Q.eff = Q.eff, par.clu = par.clu,
               VCOV = VCOV, est.gen.eff = est.gen.eff, parallel = parallel,
               cluster = cluster)

# save SIM results in output location

write.table(SIM, file = paste0(folder.loc, "/", "SIM.txt"), quote = FALSE,
            sep = "\t", row.names = FALSE)

# cofactors selection

cofactors <- QTL_select(Qprof = SIM, threshold = thre.cof, window = win.cof)


if (is.null(cofactors)) { # test if cofactors have been selected
  
  message("No QTL/cofactor position detected based on the SIM profile.")
  
  return(NULL)

  

}

# 4. Multi-QTL model search - CIM
#################################

if(!silence.print){

  cat("\n")
  cat("Multi-QTL model search - CIM")
  cat("\n")
  cat("\n")

}

CIM <- mpp_CIM(mppData = mppData, Q.eff = Q.eff, par.clu = par.clu,
               VCOV = VCOV, cofactors = cofactors, window = window,
               est.gen.eff = est.gen.eff, parallel = parallel,
               cluster = cluster)


if (N.cim > 1) {

  for (i in 1:(N.cim - 1)) {

    # take the cofactors of the previous analysis

    cofactors <- QTL_select(Qprof = CIM, threshold = thre.cof,
                            window = win.cof)

    if (is.null(cofactors)) { # test if cofactors have been selected
      
      mes.text <- paste("No cofactor position detected in CIM profile nb", i)
      message(mes.text)
      
      return(NULL)

    } else {
      
      if(!silence.print){
        
        cat("\n")
        cat(paste("CIM scan", (i+1)))
        cat("\n")
        cat("\n")
        
      }

      CIM <- mpp_CIM(mppData = mppData, Q.eff = Q.eff, par.clu = par.clu,
                     VCOV = VCOV, cofactors = cofactors, window = window,
                     est.gen.eff = est.gen.eff, parallel = parallel,
                     cluster = cluster)

    }

  }

}

# save the list of cofactors

write.table(cofactors[, 1:5], file = paste0(folder.loc, "/", "cofactors.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# save CIM results

write.table(CIM, file = paste0(folder.loc, "/", "CIM.txt"), quote = FALSE,
            sep = "\t", row.names = FALSE)

# select QTL candidates

QTL <- QTL_select(Qprof = CIM, threshold = thre.QTL, window = win.QTL)

if (is.null(QTL)) { # test if QTL have been selected

  message("No QTL position detected based on the (last) CIM profile.")
  return(NULL)

}


# 5. Backward elimination
#########################

if (backward){

  if(!silence.print){

    cat("\n")
    cat("Backward elimination")
    cat("\n")
    cat("\n")

  }

  QTL <- mpp_BackElim(mppData = mppData, QTL = QTL, Q.eff = Q.eff,
                      par.clu = par.clu, VCOV = VCOV, alpha = alpha.bk)
  
  if (is.null(QTL)) { # test if QTL have been selected
    
    stop("No QTL position stayed in the model after the backward elimination.
         This is probably due to an error in the computation of the model
         in asreml() function.")
    
  }

}

# save the final list of QTLs

write.table(QTL[, 1:5], file = paste0(folder.loc, "/", "QTL.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)


# 6. R squared computation
##########################

if(!silence.print){

  cat("\n")
  cat("R squared computation")
  cat("\n")
  cat("\n")

}


R2 <- QTL_R2(mppData = mppData, QTL = QTL, Q.eff = Q.eff, par.clu = par.clu,
             VCOV = VCOV, LR.R2 = LR.R2)

# try to calcuate R2 with the linear method if it fail with LR R2
# for VCOV = "cr.err".

if(((VCOV == "cr.err") & (LR.R2)) & is.na(R2[[1]][1])){
  
  R2 <- QTL_R2(mppData = mppData, QTL = QTL, Q.eff = Q.eff, par.clu = par.clu,
               VCOV = "h.err")
  
  message <- paste("The computation of",
                   "the likelihood R squared statistics failed probably",
                   "due to some singularities. The R2 were re-calculated",
                   "using linear models.")
  
  message(message)
  
}

# save R2 results

QTL.R2 <- data.frame(QTL[, 1:5], round(R2[[3]], 2), round(R2[[4]], 2),
                     round(R2[[5]], 2), round(R2[[6]], 2),
                     stringsAsFactors = FALSE)

colnames(QTL.R2)[6:9] <- c("R2.diff", "adj.R2.diff", "R2.sg", "adj.R2.sg")

write.table(QTL.R2, file = paste0(folder.loc, "/", "QTL_R2.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# 7. QTL effects estimation
###########################

if(!silence.print){

  cat("\n")
  cat("QTL effects estimation")
  cat("\n")
  cat("\n")

}

QTL.effects <- QTL_genEffects(mppData = mppData, QTL = QTL, Q.eff = Q.eff,
                              par.clu = par.clu, VCOV = VCOV, const = const,
                              par.ref = par.ref)

# 8. CIM- and confidence interval computation
#############################################

if(CI){

  if(!silence.print){

    cat("\n")
    cat("CIM- and confidence intervals computation")
    cat("\n")
    cat("\n")

  }

  CIM.m <- mpp_CIM(mppData = mppData, Q.eff = Q.eff, par.clu = par.clu,
                   VCOV = VCOV, cofactors = cofactors, window = 1000000,
                   est.gen.eff = FALSE, parallel = parallel,
                   cluster = cluster)

  QTL.CI <- QTL_CI(QTL = QTL, Qprof = CIM.m, drop = drop)
  
  write.table(QTL.CI, file = paste0(folder.loc, "/", "QTL_CI.txt"),
              quote = FALSE, sep = "\t", row.names = FALSE)

} else { QTL.CI <- NULL}



  # 9. Results processing
  #######################

  if(!silence.print){

    cat("\n")
    cat("Results processing")
    cat("\n")
    cat("\n")

  }

  ### 9.1: general results

  gen.res <- c(dim(QTL)[1], round(R2[[1]][1], 2), round(R2[[2]][1], 2))
  names(gen.res) <- c("nb.QTL", "glb.R2", "glb.adj.R2")

  write.table(gen.res, file = paste0(folder.loc, "/", "QTL_genResults.txt"),
            quote = FALSE, sep = "\t", col.names = FALSE)

  
  ### 9.2: Plots


   main.cim <- paste("CIM", pop.name, trait.name, Q.eff, VCOV)
   main.Qeff <- paste("QTL gen. effects", pop.name, trait.name, Q.eff, VCOV)

  if (Q.eff == "biall") {
  
  pdf(paste0(folder.loc, "/", "QTL_profile.pdf"), height = 10, width = 16)
  
  print(plot_QTLprof(Qprof = CIM, QTL = cofactors, type = "h", main = main.cim,
                     threshold = thre.QTL))
  
  dev.off()
  
  } else {
  
  # CIM profile
  
  pdf(paste0(folder.loc, "/", "QTL_profile.pdf"), height = 10, width = 16)
  
  
    print(plot_QTLprof(Qprof = CIM, QTL = cofactors, type = "l", main = main.cim,
                       threshold = thre.QTL))
  
  dev.off()
  
  # genetic effect plot
  
  if (est.gen.eff) {
    
    pdf(paste0(folder.loc, "/", "gen_eff.pdf"), height = 10, width = 16)
    
    print(plot_genEffects(Qprof = CIM, QTL = QTL, main = main.Qeff))
    
    dev.off()
    
  }
  
}

  ### 9.3: Report


if(CI) {QTL.info <- data.frame(QTL[, c(1, 2, 4, 5)], QTL.CI[, 4:8],
                               stringsAsFactors = FALSE)
} else {QTL.info <-  QTL[, c(1, 2, 4, 5)]}

QTL_report(out.file = paste0(folder.loc, "/", "QTL_REPORT.txt"),
           main = paste(pop.name, trait.name, Q.eff, VCOV), QTL.info = QTL.info,
           QTL.effects = QTL.effects, R2 = R2)
   

  ### 9.4: Return R object
   
   
   results <- list(n.QTL = dim(QTL)[1], cofactors = cofactors[, 1:5],
                   QTL = QTL[, 1:5], R2 = R2, QTL.effects = QTL.effects,
                   QTL.CI = QTL.CI)
   
   return(results)

}