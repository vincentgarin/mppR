###################
# plot_genEffects #
###################

#' plot of QTL genetic effects
#' 
#' Visualisation of the genome-wide significance of the QTL effect per cross or
#' per parents. It uses the results of a SIM or CIM profile from a
#' cross-specific, parental and ancestral model. candidate QTL positions passed
#' to the argument \code{QTL} are also drawn on the graph.
#' 
#' Plot the significance of the cross-specific or parental components of the QTL
#' effects. White colour means no significant effect. For a cross-specific
#' QTL profile: Red color means that the allele coming from parent A(1)
#' increases the phenotypic value and parent B(2) decreases it
#' and blue that parent A(1) decreases the trait and parent B(2) increases it.
#' 
#' For a parental or an ancestral QTL profile: Red colour means a significant
#' negative difference with respect to the first parent (ancestral group)
#' taken as reference, Blue represent a significant positive difference.
#' The colour intensity increase with the significance of the effect.
#' There are five colour intensities according to the p-value of the QTL
#' effect: 0.05<p-val<0.01; 0.01<p-val<0.001; 0.001<p-val<0.0001;
#' 0.0001<p-val<0.00001 and p-val< 0.00001.
#' 
#' The dashed lines indicate the position that have been introduce in the
#' argument \code{QTL}.
#' 
#' 
#' @param Qprof Object of class \code{QTLprof} returned by the function
#' \code{\link{mpp_SIM}} or \code{\link{mpp_CIM}} with argument
#' \code{est.gen.eff = TRUE}.
#' 
#' @param QTL Optional argument. Object of class \code{QTLlist} representing a
#' list of selected position obtained with the function \code{\link{QTL_select}}
#' or vector of \code{character} marker or inbetween marker positions names.
#' These positions will be plotted on the graph. Default = NULL.
#' 
#' @param main Title of the graph. Default =  "QTL genetic effects plot".
#' 
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{mpp_CIM}}, \code{\link{mpp_SIM}}, \code{\link{QTL_select}}
#' 
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData, Q.eff = "cr", est.gen.eff = TRUE)
#' QTL <- QTL_select(SIM)
#' 
#' plot_genEffects(Qprof = SIM) # without QTL positions
#' 
#' plot_genEffects(Qprof = SIM, QTL = QTL) # with QTL positions
#' 
#' @export
#' 


plot_genEffects <- function(Qprof, QTL = NULL,
                            main = "QTL genetic effects plot")
{
  
  # 1. check data format
  ######################
  
  stopifnot(inherits(Qprof, "QTLprof"))
  
  
  n.eff <- dim(Qprof)[2] - 5
  
  if(n.eff == 0) {
    
    stop("The Qprof object does not contain any QTL p-value information.
         It was probably not obtained using est.gen.eff = TRUE")
    
  }
  
  # 2. elements for the plot
  ##########################
  
  ### 2.1 list of QTL
  
  if(!is.null(QTL)){
    
    stopifnot(inherits(QTL, "QTLlist"))
    pos.Q <- QTL[, c(2, 4)]
    
  }
  
  
  ### 2.2 colour code from -5 red to 5 blue
  
  z <-  c(apply(X = Qprof[, 6:dim(Qprof)[2]], MARGIN = c(1, 2),
                FUN = color.code))
  
  ### 2.3 cross or parent indicator and chromosome
  
  y <- factor(rep(1:n.eff, each = dim(Qprof)[1]))
  
  # chr <- factor(rep(Qprof$chr, n.eff))
  
  chr <- rep(Qprof$chr, n.eff)
  
  ### 2.4 genetic map positions in cM with width between two positions.
  
  x <- rep(Qprof$pos.cM, n.eff)
  
  w <- tapply(X = Qprof$pos.cM, INDEX = Qprof$chr, FUN = function(x) c(diff(x), 1))
  w <- unlist(w)
  w <- rep(w, n.eff)
  
  pos.cM <- Qprof$pos.cM
  
  ### 2.5 legend for the y-axis (cross or parents)
  
  y.names <- colnames(Qprof)[6:dim(Qprof)[2]]
  
  ### 2.6 gather data for the plot
  
  data <- data.frame(x, y, z, chr, w)
  
  
  # 3. plot
  #########
  
  if(is.null(QTL)){ # no QTL position given
    
    pl <- ggplot(data, aes(x, y, z = z))
    pl + geom_tile(aes(fill = z, width = w)) +
      facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
      scale_fill_gradient2(limits = c(-5, 5), low = "red", mid = "white",
                           high = "blue") +
      theme_bw() + xlab("position [cM]") +
      scale_y_discrete(labels = y.names) + ggtitle(main)
    
  } else { # QTL position given
    
    pl <- ggplot(data, aes(x, y, z = z))
    pl + geom_tile(aes(fill = z, width = w)) +
      facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
      scale_fill_gradient2(limits = c(-5, 5), low = "red", mid = "white",
                           high = "blue") +
      geom_vline(aes(xintercept = pos.cM), pos.Q, linetype = "longdash") +
      theme_bw() + xlab("position [cM]") +
      scale_y_discrete(labels = y.names) + ggtitle(main)
    
  }
  
}