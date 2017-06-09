##################
# parent_cluster #
##################

#' Clustering of parental lines
#' 
#' Local clustering of the parental lines obtained using the R package
#' clushaplo (Leroux et al. 2014).
#' 
#' This function is a wrapper of \code{clusthaplo} package functions.
#' It can only be run after package Clusthaplo has been installed.
#' It can be found there:
#' \url{https://cran.r-project.org/src/contrib/Archive/clusthaplo/}. A
#' visualisation of ancestral haplotype blocks can be obtained setting
#' \code{plot = TRUE}. The plots will be saved at the location specified
#' in \code{plot.loc}. In order to cluster parental lines using hidden Markov
#' models (\code{clustering.method = "hmm"}), the user needs to use an R version
#' where package \code{clusthaplo} and \code{RHmm} can be used simultaneously.
#' R 2.14.0 is a possibility.
#' 
#' @param haplo.map Three columns \code{data.frame} with: 1)
#' marker identifiers; 2) chromosomes; and 3) \code{numeric} position in
#' centi-Morgan. This argument represents the genetic map of the parents markers
#' used for the clustering procedure. It can contain as many markers as possible
#' to increased the available information for clustering.
#' 
#' @param consensus.map \code{Data.frame} with the same format as
#' \code{haplo.map}. This is the map used for the QTL analysis.
#' 
#' @param marker.data \code{Character} marker \code{matrix} representing the
#' parent scores for the markers present in \code{haplo.map}. \strong{the
#' rownames of \code{marker.data} and the marker list of \code{haplo.map} must
#' be strictly equivalents.}
#' 
#' @param na.strings \code{Character} expression for marker scores missing value
#' (".", "NA", "-"). Default = NA.
#' 
#' @param w1 The w1 weight function in the Li&Jyang similarity score.
#' Possible values are "kernel.const", "kernel.exp", "kernel.gauss",
#' "kernel.unif", "kernel.laplace" or "kernel.null". Default = "kernel.exp".
#' 
#' @param w2 The w2 weight function in the Li&Jyang similarity score.
#' Possible values are "kernel.const", "kernel.exp", "kernel.gauss",
#' "kernel.unif", "kernel.laplace" or "kernel.null". default = "kernel.unif".
#' 
#' @param step.size \code{Numeric} value. If the length between two marker
#' position of consensus.map is bigger than step.size then clusthaplo will
#' caculate a clustering results of the parental genotype at step.size
#' interval(s). The default value is set to 10000. This means that no ancestral
#' class will be inferred in between position of the consensus map.
#' 
#' @param window \code{Numeric} value for the size of the window used for
#' clustering in centi-Morgan. The clustering procedure is done for the position
#' that is in the centre of the window taking marker scores within the window
#' into consideration.
#' 
#' @param clustering.method One of "threshold", "hmm". default value =
#' "threshold".
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
#' default value = 3.
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
#' @return Return:
#' 
#' \code{List} with the following objects
#' 
#' \item{par.clu }{\code{Integer matrix} with rows repersenting marker or
#' in between positions and column corresponding to the parents. At a single
#' marker position, parents with the same value were clustered in the same
#' ancestral group.}
#' 
#' \item{av.cl }{Average number of cluster. }
#' 
#' Plot of the ancestral haplotypes are saved in a folder at the location
#' specified in \code{plot.loc} if \code{plot = TRUE}.
#' 
#' @author Vincent Garin
#' 
#' @seealso \code{\link{parent_clusterCheck}}
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
#' \dontrun{
#' 
#' library(clusthaplo)
#' 
#' data(USNAM_geno)
#' data(USNAM_mppData)
#' data(USNAM_map)
#' 
#' # Parents genotypes and map
#' geno.par <- USNAM_geno[1:6, ]
#' geno.par <- t(geno.par)
#' map.par <- USNAM_map
#' 
#' # Map obtained after mpp.data object formation and potential
#' # inference of position in between existing markers
#' map <- USNAM_mppData$map[, -3]
#' 
#' 
#' par.clu <- parent_cluster(haplo.map = map.par, consensus.map = map,
#' marker.data = geno.par, clustering.method = "threshold",
#' step.size = 1000, window = 25, plot = TRUE, plot.loc = getwd())
#' 
#' 
#' par.clu$av.cl # Average number of inferred ancestral cluster along the genome
#' 
#' # Collect the parent clustering results for the analysis
#' 
#' par.clu <- par.clu$par.clu 
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData, Q.eff = "anc", par.clu = par.clu)
#' plot_QTLprof(SIM)
#' 
#' }
#' 
#' 
#' @export
#' 


parent_cluster <- function(haplo.map, consensus.map, marker.data,
                           na.strings = NA, w1 = "kernel.exp",
                           w2 = "kernel.unif", step.size = 10000,
                           window, K = 10, clustering.method = "threshold",
                           simulation.type = "equi", simulation.Ng = 50, 
                           simulation.Nrep = 3, threshold.quantile = 95,
                           plot = TRUE, plot.loc = getwd()) {
  
  # 1. check data format
  ######################
  
  if(!((exists("clusthaplo")) && (is.function(clusthaplo)))){
    
    stop("To use this function, you must load the package clusthaplo")
    
  }
  
  stopifnot(identical(haplo.map[, 1], rownames(marker.data)))
  
  # test if some makers are at the same position in the consensus map
  
  if ( sum((diff(consensus.map[, 3]) == 0) * 1) > 0) {
    
    stop("Some marker of the consensus.map are at the same position")
    
  }
  
  # keep marker, parents names, number of parents and chr length for later
  
  mk.names <- consensus.map[, 1]
  par.names <- colnames(marker.data)
  n.par <- length(par.names)
  
  
  # 2. transformation of the maps into MapMaker format
  ####################################################
  
  haplo.map <- toMapMaker(haplo.map)
  consensus.map <- toMapMaker(consensus.map)
  
  l.chr <- unlist(lapply(X = haplo.map, FUN = function(x) max(x[, 4])))
  
  # get the number of chromosome and the number of position per
  # chromosome for after checks.
  
  n.chr <- length(haplo.map)
  chr.ind <- paste("chr", 1:n.chr, sep = "")
  
  nb.pos <- unlist(lapply(X = consensus.map, FUN = function(x) dim(x)[1]))
  
  # 3. Clustering procedure
  #########################
  
  # Initialize the clustering procedure
  
  clust <- clusthaplo(haplotypes.map = haplo.map,
                      mcqtl.consensus.map = consensus.map, 
                      marker.data = marker.data, na.strings = na.strings,
                      discard.unknown.markers = TRUE)
  
  # change the clusthaplo configuration with the desired settings
  
  
  clust$config(w1 = w1, w2 = w2, step.size = step.size, window.length = window, 
               na.replace = NA, scoring.method = "kinship",
               clustering.method = clustering.method, 
               kinship.threshold = K, simulation.type = simulation.type,
               simulation.Ng = simulation.Ng,
               simulation.Np = dim(marker.data)[2],
               simulation.Nrep = simulation.Nrep, 
               threshold.quantile = threshold.quantile)
  
  # prepare folder path to save the plots if user want it
  
  if(plot){
    
    end.char <- substr(plot.loc, nchar(plot.loc), nchar(plot.loc))
    
    if(end.char == "/"){
      
      folder.loc <- paste0(plot.loc, "par_clu_plots")
      
    } else {
      
      folder.loc <- paste0(plot.loc, "/par_clu_plots")
      
    }
    
    dir.create(folder.loc)
    
  }
  
  tc <- c()
  
  for (i in 1:n.chr) {
    
    # produce the transitive closure
    
    clust$select.chromosome(chr.ind[i])
    clust$train()
    tc.i <- clust$transitive.closure(clust$pairwise.similarities())
    
    
    
    # check if the tc has the same length as the number of positions on the
    # chromsome
    
    if (dim(tc.i)[1] == nb.pos[i]) {
      
      tc <- rbind(tc, tc.i)
      
      if(plot){
        
        file <- paste0(folder.loc, paste0("/chr_", i, ".pdf"))
        
        pdf(file)
        
        print(plot(tc.i))
        
        dev.off()
        
      }
      
    } else {
      
      
      stop(paste("clusthaplo produce results for less position than the number",
                 "of the consensus map. Some marker of the consensus map have",
                 "potentially the same position."))
      
    }
    
    
  }
  
  rownames(tc) <- mk.names
  colnames(tc) <- par.names
  
  # 4. compute the average number of cluster
  ##########################################
  
  nb.cl <- apply(X = tc, MARGIN = 1, FUN = function(x) length(unique(x)))
  
  av.cl <- mean(nb.cl)
  
  return(list(par.clu = tc, av.cl = av.cl))
  
}