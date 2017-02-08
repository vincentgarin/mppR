##########
# mpp_CV #
##########

#' MPP cross-validation
#' 
#' Evaluation of MPP QTL detection procedure by cross-validation (CV). For
#' details on the MPP QTL detection models see \code{\link{mpp_SIM}}
#' documentation. The CV scheme is adapted from Utz et al. (2000) to the MPP
#' context.
#' 
#' A single CV run works like that:
#' 
#' \enumerate{
#' 
#' \item{Generation of a k-fold partition of the data. The partition is done
#' within crosses. Each cross is divided into k subsets. Then for the kth
#' repetition, the kth subset is used as validation set, the rest goes into the
#' training set.}
#' 
#' \item{For the kth repetition, utilisation of the training set for cofactor
#' selection and multi-QTL model search (\code{\link{mpp_SIM}} and
#' \code{\link{mpp_CIM}}). If \code{backward = TRUE}, the final list of QTLs is
#' tested simultaneously using a backward elimination
#' (\code{\link{mpp_BackElim}}).}
#' 
#' \item{Use the list of detected QTLs in the training sets to calculate
#' the proportion of genetic variance explained by all detected QTLs in the
#' training set (p.ts = R2.ts/heritability). Global p.ts (all QTLs) are returned
#' adjusted and unadjusted.
#' 
#' For each single QTL effect, difference and single R squared are also
#' calculated. Difference R squared are computed by doing the difference between
#' a model with all QTLs and a model without the ith position. Single QTL R
#' squared are obtained from a model including only a single QTL position.
#' These partial R squared are only returned unadjusted. For details about R
#' squared computation and adjustment look at \code{\link{QTL_R2}}. }
#' 
#' \item{Use the estimates of the QTL effects in the training set (B.ts) to
#' predict the phenotypic values of the validation set. y.pred.vs = X.vs*B.ts.
#' Computes the predicted R squared  in the validation set using the squared
#' Pearson correlation coefficient between the real values (y.vs) and the
#' predicted values (y.pred.vs). R2.vs = cor(y.ts,y.pred.ts)^2.
#' 
#' The predicted R squared can be computed per cross and then averaged
#' (\code{within.cross = TRUE}) or at the population level. The proportion
#' of predicted genetic variance in the validation set is obtained dividing
#' R2.vs by the heritability. Predicted R squared are only returned unadjusted
#' because we did not find so far any satisfactory solution with respect to this
#' question. 
#' 
#' Partial QTL predicted unadjusted R squared difference and single R squared
#' are also calculated.
#' 
#'   }
#' 
#' }
#' 
#' \strong{WARNING!(1)} The computation of \code{mpp_CV()} function using mixed
#' models (all models with \code{VCOV} different than \code{"h.err"}) is
#' technically possible but can be irrealistic in practice due to a reduced
#' computer power. Since a mixed model is computed at each single position it
#' can take a lot of time. From our estimation it can take between 20 to 50
#' times more time than for linear models.
#' 
#' \strong{WARNING!(2)} The estimation of the random pedigree models
#' (\code{VCOV = "pedigree" and "ped_cr.err"}) can be unstable. Sometimes the
#' \code{asreml()} function fails to produce a results and returns the following
#' message: \strong{\code{GIV matrix not positive definite: Singular pivots}}.
#' So far we were not able to identify the reason of this problem and to
#' reproduce this error because it seems to happen randomly. The consequence of
#' this is that part of the results of the CV procedure can not be produced.
#' These will be replaced by NA values. According to the tests we run this
#' can affect until 10 percent of the results.
#' 
#' @param pop.name \code{Character} name of the studied population.
#' Default = "MPP_CV".
#' 
#' @param trait.name \code{Character} name of the studied trait.
#' Default = "trait1".
#'
#' @param mppData An object of class \code{mppData}. See
#' \code{\link{mppData_form}} for details.
#' 
#' @param her \code{Numeric} value between 0 and 1 representing the heritability
#' of the trait. \strong{By default, the heritability is set to 1
#' (\code{her = 1}). This means that the results represent the proportion of
#' phenotypic variance explained (predicted) in the training (validation) sets.}
#' 
#' @param Rep \code{Numeric} value representing the number of repetition of the
#' k-fold procedure. Default = 10.
#' 
#' @param k \code{Numeric} value representing the number of folds for the within
#' cross partition of the population. Default = 5.
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
#' @param within.cross \code{Logical} value indicating if the predicted
#' R squared must be computed within cross. In that case, predicted R squared
#' are computed within cross and the average is returned. Default = TRUE.
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
#' \code{mpp_CV()} function must be silenced. It will not
#' affect the printing of the other functions called by \code{mpp_CV()},
#' especially the printing of \code{asreml()}. Default = FALSE.
#'
#' @param output.loc Path where a folder will be created to save the results.
#' By default the function uses the current working directory.
#' 
#'    
#' @return 
#' 
#' \code{List} containing the following results items:
#' 
#' \item{QTL}{\code{Data.frame} with as row the QTL position detected at least
#' 1 time during the entire CV process and as column:
#' 
#' \enumerate{
#' 
#' \item{QTL position information (chromosome, position [cM], etc.).}
#' \item{Number of time the position has been detected (N).}
#' 
#' \item{Average  proportion of explained genetic variance in the training set
#' difference R squared (av.pts.d).}
#' \item{Average  proportion of predicted genetic variance in the validation set
#' difference R squared (av.pvs.d).}
#' \item{Relative bias explained genetic variance difference R squared (bias.d).}
#' 
#' \item{Average  proportion of explained genetic variance in the training set
#' single QTL R squared (av.pts.s).}
#' \item{Average  proportion of predicted genetic variance in the validation set
#' single QTL R squared (av.pvs.s).}
#' \item{Relative bias explained genetic variance single QTL R squared (bias.s).}
#'
#' }
#' 
#' } 
#' 
#' \item{N.QTL}{Table with the number of QTL detected for each (Rep*k)
#' repetition.}
#' 
#' \item{p.ts}{Table with the proportion of explained genetic variance from
#' unadjusted R squared in the training set for each (Rep*k) repetition.}
#' 
#' \item{p.vs}{Table with the proportion of predicted genetic variance from
#' unadjusted R squared in the validation set for each (Rep*k) repetition.}
#' 
#' \item{bias}{Explained genetic variance bias between training and validation
#' set from unadjusted R squared 1-(p.vs/p.ts).}
#' 
#' \item{p.ts.adj}{Table with the proportion of explained genetic variance from
#' adjusted R squared in the training set for each (Rep*k) repetition.}
#' 
#' \item{QTL.profiles}{Final QTL profiles of each (Rep*k) repetition.}
#' 
#' The above mentioned elements return as R object are also saved  as text
#' files at the specified output location (\code{output.loc}). A transparency
#' plot of the CV results (plot.pdf) obtained with the function
#' \code{\link{plot_CV}} is also saved.
#' 
#' 
#' @author Vincent Garin
#' 
#' @references
#' 
#' Utz, H. F., Melchinger, A. E., & Schon, C. C. (2000). Bias and sampling error
#' of the estimated proportion of genotypic variance explained by quantitative
#' trait loci determined from experimental data in maize using cross validation
#' and validation with independent samples. Genetics, 154(4), 1839-1849.
#' 
#' @seealso
#' 
#' \code{\link{mppData_form}},
#' \code{\link{mpp_BackElim}},
#' \code{\link{mpp_CIM}},
#' \code{\link{mpp_perm}},
#' \code{\link{mpp_SIM}},
#' \code{\link{parent_cluster}},
#' \code{\link{plot_CV}},
#' \code{\link{QTL_CI}},
#' \code{\link{QTL_R2}},
#' \code{\link{QTL_pred_R2}},
#' \code{\link{USNAM_parClu}}
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(USNAM_mppData)
#' 
#' my.loc <- "C:/..."
#' 
#' CV <- mpp_CV(pop.name = "USNAM", trait.name = "ULA", mppData = USNAM_mppData,
#' her = .5, output.loc = my.loc)
#' 
#' plot_CV(CV.res = CV)
#' 
#' # Using parallelization
#' 
#' library(parallel)
#' n.cores <- detectCores()
#' cluster <- makeCluster(n.cores-1)
#' 
#' CV <- mpp_CV(pop.name = "USNAM", trait.name = "ULA", mppData = USNAM_mppData,
#'              her = .5, output.loc = my.loc, parallel = TRUE, cluster = cluster)
#' 
#' plot_CV(CV.res = CV)
#' 
#' }
#' 
#' @export
#'


mpp_CV <- function(pop.name = "MPP_CV", trait.name = "trait1",
                   mppData, her = 1, Rep = 10, k = 5, Q.eff = "cr",
                   par.clu = NULL, VCOV = "h.err", thre.cof = 3, win.cof = 20,
                   N.cim = 1, window = 20, thre.QTL = 3, win.QTL = 20,
                   backward = TRUE, alpha.bk = 0.05, within.cross = TRUE,
                   parallel = FALSE, cluster = NULL, silence.print = FALSE,
                   output.loc = getwd()) {
  
  # 1. Check the validity of the parameters that have been introduced
  ###################################################################
  
  check.mpp.cv(mppData = mppData, Q.eff = Q.eff, VCOV = VCOV,
               par.clu = par.clu, parallel = parallel, cluster = cluster,
               output.loc = output.loc)
  
  # 2. Create a directory to store the results
  ############################################
  
  # create a directory to store the results of the QTL analysis
  
  end.char <- substr(output.loc, nchar(output.loc), nchar(output.loc))
  
  if(end.char == "/"){
    
    folder.loc <- paste0(output.loc, paste("CV",pop.name, trait.name, Q.eff,
                                           VCOV, sep = "_"))
    
  } else {
    
    folder.loc <- paste0(output.loc, "/", paste("CV",pop.name, trait.name,
                                                Q.eff, VCOV, sep = "_"))
    
  }
  
  dir.create(folder.loc)
  
  # 3. Create space to store the results
  ######################################
  
  N.QTL <- p.ts <- p.ts.adj <- p.vs <- bias <- matrix(0, k, Rep)
  
  QTL.positions <- rep(0, dim(mppData$map)[1])
  N.Qeff.est.ts <- rep(0, dim(mppData$map)[1])
  N.Qeff.est.vs <- rep(0, dim(mppData$map)[1])
  
  QTL.pts.s <- QTL.pvs.s <- rep(0, dim(mppData$map)[1])
  
  QTL.pts.d <- QTL.pvs.d <- rep(0, dim(mppData$map)[1])
  
  profiles <- c()
  
  # keep the marker and in between position full list
  
  mk.list <- mppData$map[, 1]
  
  
  # 4. start to loop from 1 to r replicates
  ########################################
  
  for (i in 1:Rep) {
    
    if(!silence.print){
      
      cat(paste("CV repetition", i))
      cat("\n")
      
    }
    
    
    ### 4.1 generate a CV partition
    
    folds <- CV_partition(cross.ind = mppData$cross.ind, k = k)
    
    ### 4.2 iterate through the CV partition
    
    for (j in 1:k) {
      
      if(!silence.print){
        
        cat(paste("fold", j))
        cat("\n")
        
      }
      
      # training set
      
      mppData.ts <- mppData_subset(mppData = mppData,
                                   gen.list = folds[[j]]$train.set)
      
      # validation set
      
      mppData.vs <- mppData_subset(mppData = mppData,
                                   gen.list = folds[[j]]$val.set)
      
      prob.prog <- FALSE # indicator variable for a programmation problem
      
      # 4.2.1 cofactors selection
      
      SIM <- mpp_SIM(mppData = mppData.ts, Q.eff = Q.eff, par.clu = par.clu,
                     VCOV = VCOV, parallel = parallel, cluster = cluster)
      
      if(sum(SIM$log10pval) == 0){prob.prog <- TRUE }
      
      cofactors <- QTL_select(Qprof = SIM, threshold = thre.cof,
                              window = win.cof, silence.print = TRUE)
      
      
      # 4.2.2 multi-QTL model search
      
      # test if some cofactors have been selected
      
      if (!is.null(cofactors)) {
        
        # there are some cofactors
        
        CIM <- mpp_CIM(mppData = mppData.ts, Q.eff = Q.eff,
                       par.clu = par.clu, VCOV = VCOV, cofactors = cofactors,
                       window = window, parallel = parallel, cluster = cluster)
        
        if(sum(CIM$log10pval) == 0){prob.prog <- TRUE }
        
        if (N.cim > 1) {
          
          for (l in 1:(N.cim - 1)) {
            
            # take the cofactors of the previous analysis
            
            cofactors <- QTL_select(Qprof = CIM, threshold = thre.cof,
                                    window = win.cof, silence.print = TRUE)
            
            # test if some cofactors there before running next CIM
            
            if (!is.null(cofactors)) {
              
              CIM <- mpp_CIM(mppData = mppData.ts, Q.eff = Q.eff,
                             par.clu = par.clu, VCOV = VCOV,
                             cofactors = cofactors,
                             window = window, parallel = parallel,
                             cluster = cluster)
              
              if(sum(CIM$log10pval) == 0){prob.prog <- TRUE }
              
              # else leave the loop
              
            } else { break }
            
          }
          
        }
        
        
        
        #### end multi QTL search
        
        # 4.2.3 QTL selection
        
        QTL <- QTL_select(Qprof = CIM, threshold = thre.QTL, window = win.QTL,
                          silence.print = TRUE)
        
        # 4.2.4 Optional backward elimination
        
        if(!is.null(QTL) & backward){
          
          QTL.back <- mpp_BackElim(mppData = mppData.ts, QTL = QTL, Q.eff = Q.eff,
                                   par.clu = par.clu, VCOV = VCOV, alpha = alpha.bk)
          
          if(is.null(QTL.back)){ # If there was QTL position and backward return
            # no QTL it is (probably) due to programming error.
            
            prob.prog <- TRUE
            QTL <- NULL
            
          } else {QTL <- QTL.back}
          
        }
        
        
        if (!is.null(QTL)) {
          
          
          ### 4.3 compute the CV statistics (N.QTL, p.ts. etc.)
          
          # a) N.QTL
          
          N.QTL[j, i] <- dim(QTL)[1]
          
          # b) QTL positions
          
          QTL.names <- QTL[, 1]
          
          QTL.positions <- QTL.positions + (is.element(mk.list, QTL.names) * 1)
          
          # store the profiles
          
          profiles <- cbind(profiles, CIM$log10pval)
          
          # c) p.ts (adjusted R2 / heritability)
          
          # get the R squared training set
          ################################
          
          
          R2.ts <- QTL_R2(mppData = mppData.ts, QTL = QTL, Q.eff = Q.eff,
                          par.clu = par.clu)
          
          
          # compute predicted R squared
          ##############################
          
          R2.vs <- QTL_pred_R2(mppData.ts = mppData.ts, mppData.vs = mppData.vs,
                               Q.eff = Q.eff, par.clu = par.clu, VCOV = VCOV,
                               QTL = QTL, within.cross = within.cross)
          
          
          # global results
          ################
          
          # non adjusted
          
          p.ts[j, i] <- R2.ts$glb.R2/her
          p.vs[j, i] <- R2.vs$glb.R2/her
          bias[j, i] <- round(1 - (p.vs[j, i]/p.ts[j, i]), 2)
          
          
          # adjusted
          
          p.ts.adj[j, i] <- R2.ts$glb.adj.R2/her
          
          
          # store partial R squared
          #########################
          
          # count the number of time partial QTL R2 could be estimated 
          
          if(!is.na(R2.vs[[1]])){ # validation set
            
            N.Qeff.est.vs <- N.Qeff.est.vs + (is.element(mk.list, QTL.names) * 1)
            
          }
          
          if(!is.na(R2.ts[[1]])){ # test set
            
            N.Qeff.est.ts <- N.Qeff.est.ts + (is.element(mk.list, QTL.names) * 1)
            
          }
          
          
          # general function to store individual QTL results
          
          update.res <- function(input, output, Q.names, mk.list){
            
            R.ij <- rep(0, length(mk.list))
            R.ij[match(Q.names, mk.list)] <- input
            rowSums(data.frame(output, R.ij), na.rm = TRUE)
            
          }
          
          # single QTL R2 values
          
          QTL.pts.s <- update.res(input = R2.ts$part.R2.sg/her,
                                  output = QTL.pts.s, Q.names = QTL.names,
                                  mk.list = mk.list)
          
          QTL.pvs.s <- update.res(input = R2.vs$part.R2.sg/her,
                                  output = QTL.pvs.s, Q.names = QTL.names,
                                  mk.list = mk.list)
          
          # difference QTL R2 values
          
          QTL.pts.d <- update.res(input = R2.ts$part.R2.diff/her,
                                  output = QTL.pts.d, Q.names = QTL.names,
                                  mk.list = mk.list)
          
          QTL.pvs.d <- update.res(input = R2.vs$part.R2.diff/her,
                                  output = QTL.pvs.d, Q.names = QTL.names,
                                  mk.list = mk.list)
          
        } else {
          
          
          if(prob.prog){
            
            N.QTL[j, i] <- p.ts[j, i] <- p.vs[j, i] <- bias[j, i] <- NA
            p.ts.adj[j, i] <- NA
            
          } else { # only no QTL significant enough.
            
            # N.QTL p.ts and p.vs stay at 0 as in the initialisation. Only need to
            # put the bias at NA.
            
            bias[j, i] <- NA
            
            profiles <- cbind(profiles, CIM$log10pval)
            
          }
          
        }
        
      } else {
        
        
        if(prob.prog){
          
          N.QTL[j, i] <- p.ts[j, i] <- p.vs[j, i] <- bias[j, i] <- NA
          p.ts.adj[j, i] <- NA
          
        } else { # Only no QTL significant enough.
          
          # N.QTL p.ts and p.vs stay at 0 as in the initialisation. Only need to
          # put the bias at NA.
          
          bias[j, i] <- NA
          
          profiles <- cbind(profiles, SIM$log10pval)
          
        }
        
      }
      
    }  # end jth fold loop
    
    
  }  # end ith repetition loop
  
  
  # 5. save the results in the specified folder
  #############################################
  
  
  save.CV.res <- function(res, col, row, out, name){
    colnames(res) <- col; rownames(res) <- row
    write.table(res, paste0(out, "/", name, ".txt"), quote = FALSE, sep = "\t")
    
  }
  
  
  row.n <- paste0("k_", 1:k)
  col.n <- paste0("Rep_", 1:Rep)
  
  save.CV.res(res = N.QTL, col = col.n, row = row.n, out = folder.loc,
              name = "N_QTL")
  
  save.CV.res(res = p.ts, col = col.n, row = row.n, out = folder.loc,
              name = "p_ts")
  
  save.CV.res(res = p.ts.adj, col = col.n, row = row.n, out = folder.loc,
              name = "p_ts_adj")
  
  save.CV.res(res = p.vs, col = col.n, row = row.n, out = folder.loc,
              name = "p_vs")
  
  save.CV.res(res = bias, col = col.n, row = row.n, out = folder.loc,
              name = "bias")
  
  
  # QTL.positions
  
  av.pts.s <- round(QTL.pts.s/N.Qeff.est.ts, 1)
  av.pvs.s <- round(QTL.pvs.s/N.Qeff.est.vs, 1)
  bias.s <- round(1 - (av.pvs.s/av.pts.s), 1)
  
  av.pts.d <- round(QTL.pts.d/N.Qeff.est.ts, 1)
  av.pvs.d <- round(QTL.pvs.d/N.Qeff.est.vs, 1)
  bias.d <- round(1 - (av.pvs.d/av.pts.d), 1)
  
  
  QTL.sum <- data.frame(mppData$map, QTL.positions, av.pts.d, av.pvs.d,
                        bias.d, av.pts.s, av.pvs.s, bias.s,
                        stringsAsFactors = FALSE)
  
  QTL.sum <- QTL.sum[QTL.positions > 0, ]
  
  colnames(QTL.sum)[5] <- "N"
  
  write.table(QTL.sum, paste0(folder.loc, "/", "QTL.txt"), quote = FALSE,
              sep = "\t", row.names = FALSE)
  
  # profiles
  
  file <- paste0(folder.loc, "/", "QTL_profiles.txt")
  
  # combine the profiles with the map
  
  profiles2 <- data.frame(mppData$map, QTL.positions, profiles,
                          stringsAsFactors = FALSE)
  
  colnames(profiles2)[5:dim(profiles2)[2]] <- c("N.QTL", paste0("QTL.prof_",
                                                                1:dim(profiles)[2]))
  
  write.table(profiles2, file, quote = FALSE, sep = "\t", row.names = FALSE)
  
  # produce profile plot
  
  pdf(paste0(folder.loc, "/", "plot.pdf"), height = 10, width = 16)
  
  plot_CV(CV.res = list(QTL.profiles = profiles2),
          main = paste("CV", trait.name, Q.eff, VCOV))
  
  dev.off()
  
  # return R object
  
  return(list(QTL = QTL.sum, N.QTL = N.QTL, p.ts = p.ts, p.vs = p.vs,
              bias = bias, p.ts.adj = p.ts.adj, QTL.profiles = profiles2))
  
  
} 
