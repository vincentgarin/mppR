##########################
# parent_cluster.mppData #
##########################

#' Parent clustering for \code{mppData} objects
#' 
#' Local clustering of the parental lines done by the R package clushaplo
#' (Leroux et al. 2014) or by providing own parent clustering data.
#' 
#' This function integrate the parent clustering information to the mppData
#' object. The parent clustering is necessary to compute the ancestral model.
#' If the parent clustering step is skipped, the ancestral model can not be
#' used but the other models (cross-specific, parental, and bi-allelic) can
#' still be computed.
#' 
#' The parent clustering can be performed using the R package
#' clusthaplo using \code{method = "clusthaplo"}. Clusthaplo can be found there:
#' \url{https://cran.r-project.org/src/contrib/Archive/clusthaplo/}. Using
#' clusthaplo, a visualisation of ancestral haplotype blocks can be obtained
#' setting \code{plot = TRUE}. The plots will be saved at the location specified
#' in \code{plot.loc}.
#' 
#' An alternative (\code{method = "given"}), is to provide your own parent
#' clustering information via the argument \code{par.clu}.
#' 
#' @param mppData  An object of class \code{mppData}. the \code{mppData} must
#' have been processed using: \code{\link{create.mppData}},
#' \code{\link{QC.mppData}}, \code{\link{IBS.mppData}},
#' and \code{\link{IBD.mppData}}.
#' 
#' @param method \code{Character} expression. If (\code{method = "clusthaplo"}),
#' the clustering is done using the R package clusthaplo.
#' If (\code{method = "given"}), the user must provide the parent clustering
#' information using \code{par.clu}. Default = NULL.
#' 
#' @param par.clu Optional argument if (\code{method = "given"}).
#' \code{Interger matrix} representing the results of a
#' parents genotypes clustering. The columns represent the parental lines and
#' the rows the markers. The columns names must be the same as the parents
#' list of the mppData object. The rownames must be the same as the map marker
#' list of the mppData object. At a particular position, parents with the same
#' value are assumed to inherit from the same ancestor. for more details,
#' see \code{\link{par_clu}}. Default = NULL.
#' 
#' @param w1 The w1 weight function in the Li&Jyang similarity score.
#' Possible values are "kernel.const", "kernel.exp", "kernel.gauss",
#' "kernel.unif", "kernel.laplace" or "kernel.null". Default = "kernel.exp".
#' 
#' @param w2 The w2 weight function in the Li&Jyang similarity score.
#' Possible values are "kernel.const", "kernel.exp", "kernel.gauss",
#' "kernel.unif", "kernel.laplace" or "kernel.null". Default = "kernel.unif".
#' 
#' @param window \code{Numeric} value for the size of the window used for
#' clustering in centi-Morgan. The clustering procedure is done for the position
#' that is in the centre of the window taking marker scores within the window
#' into consideration.
#' 
#' @param K A positive integer representing the number of markers in a window
#' below which the kinship data will be used. Default = 10.
#' 
#' @param simulation.type The type of simulation used for the training.
#' One of "equi" or "mosaic". Default = "equi".
#' 
#' @param simulation.Ng The number of intermediary generations to simulate
#' for the training (only relevant for "mosaic"). Default = 50.
#' 
#' @param simulation.Nrep The number of replicates to simulate for the training.
#' Default = 3.
#' 
#' @param threshold.quantile The quantile to use to select the threshold
#' automatically. It must be a plain integer xx with 80 <= xx <= 100.
#' Default = 95.
#' 
#' @param plot \code{Logical} value indicating if the plot of the clustering
#' results must be saved at the location specified in argument \code{plot.loc}.
#' Default = TRUE.
#' 
#' @param plot.loc Path where a folder will be created to save the plot of
#' the clustering results. By default the function uses the current working
#' directory.
#' 
#' 
#' 
#' @return
#' 
#' an increased \code{mppData} object containing the the same elements
#' as the \code{mppData} object provided as argument and the
#' following new elements:
#' 
#' \item{par.clu}{\code{Integer matrix} with rows repersenting markers and
#' columns corresponding to the parents. At a single marker position, parents
#' with the same value were clustered in the same ancestral group.}
#' 
#' \item{n.anc}{Average number of ancestral clusters along the genome.}
#' 
#' \item{mono.anc}{Positions for which the ancestral clustering was monomorphic.}
#' 
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{create.mppData}}, \code{\link{QC.mppData}},
#' \code{\link{IBS.mppData}}, \code{\link{IBD.mppData}}
#' 
#' @references
#' 
#' Leroux, D., Rahmani, A., Jasson, S., Ventelon, M., Louis, F., Moreau, L.,
#' & Mangin, B. (2014). Clusthaplo: a plug-in for MCQTL to enhance QTL detection
#' using ancestral alleles in multi-cross design. Theoretical and Applied
#' Genetics, 127(4), 921-933. 
#' 
#' @examples
#' 
#' data(mppData_init)
#' data(par_clu)
#' 
#' mppData <- QC.mppData(mppData_init)
#' mppData <- IBS.mppData(mppData = mppData)
#' 
#' mppData <- IBD.mppData(mppData = mppData, type = 'RIL',
#'                        type.mating = 'selfing')
#'                        
#' mppData <- parent_cluster.mppData(mppData = mppData, method = "given",
#'                                   par.clu  = par_clu)                         
#'                        
#' \dontrun{
#'                                                 
#' library(clusthaplo)
#'                         
#' mppData <- parent_cluster.mppData(mppData = mppData, method = "clusthaplo",
#'                                   window = 25, K = 10, plot = FALSE)                        
#' 
#' }
#' 
#' @export

parent_cluster.mppData <- function(mppData, method = NULL, par.clu = NULL,
                                   w1 = "kernel.exp", w2 = "kernel.unif",
                                   window, K = 10, simulation.type = "equi",
                                   simulation.Ng = 50,  simulation.Nrep = 3,
                                   threshold.quantile = 95, plot = TRUE,
                                   plot.loc = getwd()){
  
  # check the format of the data
  #################################
  
  if(!is_mppData(mppData)){
    
    stop("'mppData' must be of class ", dQuote("mppData"))
    
  }
  
  # test if correct step in the mppData processing
  
  if(!(mppData$status %in% c('IBD', 'complete'))){
    
    stop("you have to process 'mppData' in a strict order: ",
         "create.mppData, QC.mppData, IBS.mppData, IBD.mppData, ",
         "parent_cluster.mppData. You can only use parent_cluster.mppData ",
         "after create.mppData, QC.mppData, IBS.mppData, and IBD.mppData")
    
  }
  
  # check method
  
  if (is.null(method)){
    
    stop("'method' is not provided")
    
  }
  
  if(!(method %in% c("clusthaplo", "given"))){
    
    stop("'method' must be ", dQuote("clusthaplo"), ' or ', dQuote("given"))
    
  }
  
  if (method == "clusthaplo") {
    
    # Test if clusthaplo is available
    #################################
    
    test <- requireNamespace(package = 'clusthaplo', quietly = TRUE)
    
    if(!test){
      
      stop("the clusthaplo library is not available")
      
    }
    
    # Restore the necessary objects from the mppData object
    ########################################################
    
    haplo.map <- mppData$haplo.map
    consensus.map <- mppData$map[, -3]
    map <- mppData$map
    marker.data <- t(mppData$geno.par.clu)
    parents <- mppData$parents
    
    # cluster the parents
    #####################
    
    # For the step size, determine the value of the largest chromsome
    
    chr.fact <- factor(x = map[, 2], levels = unique(map[, 2]))
    
    step.size <- max(tapply(X = map[, 4], INDEX = chr.fact, FUN = max)) + 100 
    
    p_clu <- parent_cluster(haplo.map = haplo.map, consensus.map = consensus.map,
                            marker.data = marker.data, na.strings = NA,
                            w1 = w1, w2 = w2, step.size = step.size,
                            window = window, K = K,
                            simulation.type = simulation.type,
                            simulation.Ng = simulation.Ng,
                            simulation.Nrep = simulation.Nrep,
                            threshold.quantile = threshold.quantile, plot = plot,
                            plot.loc = plot.loc)
    
    # Check the monomorphic positions
    #################################
    
    par.clu <- parent_clusterCheck(par.clu = p_clu[[1]])
    
    # put the column order as the parents
    
    p_c <- par.clu[[1]]
    p_c <- p_c[, parents]
    
    # Fill the mppData object
    ##########################
    
    mppData$par.clu <- p_c
    
    nb.cl <- apply(X = p_c, MARGIN = 1, FUN = function(x) length(unique(x)))
    
    mppData$n.anc <- mean(nb.cl)
    
    mppData$mono.anc <- par.clu[[2]]
    
    mppData$status <- 'complete'
    
    mppData$geno.off <- NULL
    
    class(mppData) <- c("mppData", "list")
    
    return(mppData)
    
    
  } else { # method = "given"
    
    if(!is.matrix(par.clu)){
      
      stop("'par.clu' argument is not a matrix")
      
    }
    
    if(!is.integer(par.clu)){
      
      stop("'par.clu' is not integer")
      
    }
    
    # list parent
    
    new_par <- colnames(par.clu)
    
    if(!all(new_par %in% mppData$parents)) {
      
      wrong.par <- new_par[!(new_par %in% mppData$parents)]
      pbpar <- paste(wrong.par, collapse = ", ")
      
      message <- sprintf(ngettext(length(wrong.par),
                                  "the following parent %s is not present in 'mppData'",
                                  "the following parents %s are not present in 'mppData'"),
                         pbpar)
      
      stop(message)
      
    }
    
    # list markers
    
    if(!identical(rownames(par.clu), mppData$map[, 1])){
      
      stop("the markers of 'par.clu' and 'mppData' are not identical")
      
    }
    
    # Check monomorphism in par.clu
    ###############################
    
    par.clu <- par.clu[, mppData$parents]
    
    par_clu <- parent_clusterCheck(par.clu = par.clu)
    
    # Calculate the number of ancestral cluster
    ###########################################
    
    nb.cl <- apply(X = par_clu[[1]], MARGIN = 1,
                   FUN = function(x) length(unique(x)))
    
    av.cl <- mean(nb.cl)
    
    mppData$par.clu <- par_clu[[1]]
    
    mppData$n.anc <- av.cl
    
    mppData$mono.anc <- par_clu[[2]]
    
    mppData$status <- 'complete'
    
    mppData$geno.off <- NULL
    
    class(mppData) <- c("mppData", "list")
    
    return(mppData)
    
    
    
  }
  
  
}