##########################
# parent_cluster.mppData #
##########################

#' Parent clustering for \code{mppData} objects
#' 
#' Local clustering of the parental lines done by the R package clushaplo
#' (Leroux et al. 2014).
#' 
#' This function is a wrapper for \code{clusthaplo} R package functions.
#' It can only be run after \code{clusthaplo} has been installed.
#' Clusthaplo can be found there:
#' \url{https://cran.r-project.org/src/contrib/Archive/clusthaplo/}. A
#' visualisation of ancestral haplotype blocks can be obtained setting
#' \code{plot = TRUE}. The plots will be saved at the location specified
#' in \code{plot.loc}.
#' 
#' @param mppData  An object of class \code{mppData}. the \code{mppData} must
#' have been processed using: \code{\link{create.mppData}},
#' \code{\link{QC.mppData}}, \code{\link{IBS.mppData}},
#' and \code{\link{IBD.mppData}}.
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
#' library(clusthaplo)
#' data(mppData_init)
#' 
#' mppData <- QC.mppData(mppData_init)
#' mppData <- IBS.mppData(mppData = mppData)
#' 
#' \dontrun{
#' 
#' mppData <- IBD.mppData(mppData = mppData, type = 'RIL',
#'                        type.mating = 'selfing')
#'                        
#' mppData <- parent_cluster.mppData(mppData = mppData, window = 25, K = 10,
#'                                   plot = FALSE)                        
#' 
#' }           
#' 
#' @export

parent_cluster.mppData <- function(mppData, w1 = "kernel.exp",
                                   w2 = "kernel.unif", window, K = 10,
                                   simulation.type = "equi", simulation.Ng = 50, 
                                   simulation.Nrep = 3, threshold.quantile = 95,
                                   plot = TRUE, plot.loc = getwd()){
  
  # 1. check the format of the data
  #################################
  
  if(!is_mppData(mppData)){
    
    stop(paste('the mppData provided provided is not a mppData object.',
               'Please use function create.mppData().'))
    
  }
  
  # test if correct step in the mppData processing
  
  if(mppData$status != 'IBD'){
    
    stop(paste('You have to process the mppData objects in a strict order:',
               'create.mppData(), QC.mppData(), IBS.mppData(), IBD.mppData(),',
               'parent_cluster.mppData(). You can only use parent_cluster.mppData()',
               'after performing create.mppData(), QC.mppData(),',
               'IBS.mppData(), and IBD.mppData().'))
    
  }
  
  # the check of the other argument is done in parent_cluster function
  
  
  # 2. Restore the necessary objects from the mppData object
  ##########################################################
  
  haplo.map <- mppData$haplo.map
  consensus.map <- mppData$map[, -3]
  map <- mppData$map
  marker.data <- t(mppData$geno.par.clu)
  parents <- mppData$parents
  
  # 3. cluster the parents
  ########################
  
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
  
  # 4. Check the monomorphic positions
  ####################################
  
  par.clu <- parent_clusterCheck(par.clu = p_clu[[1]])
  
  # put the column order as the parents
  
  p_c <- par.clu[[1]]
  p_c <- p_c[, parents]
 
  # 5. fill the mppData object
  #############################
  
  mppData$par.clu <- p_c
  
  mppData$n.anc <- p_clu[[2]]
  
  mppData$mono.anc <- par.clu[[2]]
  
  mppData$status <- 'complete'
  
  mppData$geno.off <- NULL
  
  class(mppData) <- c("mppData", "list")
  
  return(mppData)
  
}