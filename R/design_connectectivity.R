########################
# design_connectivity #
########################

#'
#' Determines the connected parts of a design
#' 
#' Uses the metod of Weeks and Williams (1964) and the package igraph to
#' determine connected parts of the design.
#' 
#' @param par.per.cross Three columns \code{character matrix} specifying :
#' 1) the cross indicators ; 2) the parents 1 identifiers of the crosses;
#' 3) the parents 2 identifiers of the crosses.
#' 
#' @param plot.des \code{Logical} value indicating if connected part of the
#' design should be plotted. Default = TRUE.
#' 
#' @param output.loc Path where the plot of the design will be saved if the
#' argument is given. Default = NULL.
#' 
#' @return
#' 
#' Return a list with each element representing one connected part of the design
#' and the list of parents contained in this part.
#' 
#' @author Vincent Garin
#' 
#' @references 
#' 
#' Weeks, D. L., & Williams, D. R. (1964). A note on the determination of
#' connectedness in an N-way cross classification. Technometrics, 6(3), 319-324.
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' 
#' par.per.cross <- USNAM_mppData$par.per.cross
#' 
#' con.part <- design_connectivity(par.per.cross)
#' 
#' @export
#' 

design_connectivity <- function(par.per.cross, plot.des = TRUE,
                              output.loc = NULL){
  
  # 1. Test format of the arguments
  #################################
  
  if(!is.matrix(par.per.cross)){
    
    stop("The par.per.cross argument is not a matrix")
  }
  
  if(!is.character(par.per.cross)){
    
    stop("The par.per.cross argument is not a character matrix")
    
  }
  
  if(!is.null(output.loc)){
    
    if(!file.exists(output.loc)){
      
      stop("The path specified in the argument output.loc is not valid.")
      
    }
    
  }
  
  if((!is.null(output.loc) && !plot.des)){ plot.des <- TRUE }
  
  # 2. Determine design conectedness and plot
  ###########################################
  
  vertices <- apply(X = par.per.cross[, 2:3], MARGIN = 1, FUN = function(x) x)
  
  vertices <- c(vertices)
  
  g <- graph(vertices)
  
  if(plot.des) {plot(g)}
  
  clu <- components(g)
  
  grp <- groups(clu)
  
  return(grp)
  
  
}