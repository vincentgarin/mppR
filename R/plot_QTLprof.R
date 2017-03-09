################
# plot_QTLprof #
################

#' plot QTL profile
#' 
#' Plots the -log10(p-val) profile of a QTL analysis using the package ggplot2.
#' The user can pass a list of cofactors or QTL position to the argument
#' \code{QTL}. These positions will be drawn on the graph using dotted lines.
#' 
#' @param Qprof Object of class \code{QTLprof} returned by the function
#' \code{\link{mpp_SIM}} or \code{\link{mpp_CIM}}.
#' 
#' @param QTL Optional argument. Object of class \code{QTLlist} representing a
#' list of selected position obtained with the function \code{\link{QTL_select}}.
#' These positions will be indicated on the graph. Default = NULL.
#' 
#' @param type \code{Character} expression indicating the type of plot should be
#' drawn: "l" for lines , "h" for vertical bar. Default = "l".
#' 
#' @param main Title of the graph. Default = "QTL profile".
#' 
#' @param threshold \code{Numeric} QTL significance threshold value draw on
#' the plot. Default = 3.
#' 
#' @author Vincent Garin
#' 
#' @seealso
#' 
#' \code{\link{mpp_SIM}}, \code{\link{mpp_CIM}}, \code{\link{QTL_select}}
#' 
#' @examples
#' 
#' data(USNAM_mppData)
#' data(USNAM_mppData_bi)
#' 
#' SIM <- mpp_SIM(mppData = USNAM_mppData, Q.eff = "cr")
#' QTL <- QTL_select(SIM)
#' plot_QTLprof(Qprof = SIM, QTL = QTL)
#' 
#' SIM_biall <- mpp_SIM(mppData = USNAM_mppData_bi, Q.eff = "biall")
#' plot_QTLprof(Qprof = SIM_biall, type = "h")
#' 
#' @export
#' 


plot_QTLprof <- function(Qprof, QTL = NULL, type = "l", main = "QTL profile",
                         threshold = 3)
{
  stopifnot(inherits(Qprof, "QTLprof"))
  
  # QTL positions
  
  if(!is.null(QTL)){
    
    stopifnot(inherits(QTL, "QTLlist"))
    pos.Q <- QTL[, c(2, 4)]
    
  }
  
  # redefine data within the function to suppress R CMD check notation
  
  chr <- Qprof$chr
  pos.cM <- Qprof$pos.cM
  log10pval <- Qprof$log10pval
  Qprof <- data.frame(chr, pos.cM, log10pval)
  
  
  if(is.null(QTL)){ # no QTL info given
    
    if (type == "l") {
      
      ggplot(Qprof, aes(x = pos.cM, y = log10pval, group = chr)) + geom_line() + 
        facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
        geom_hline(yintercept = threshold, colour = "red") + theme_bw() +
        xlab("position [cM]") + ylab("-log10(p.val)") + 
        ggtitle(main) +
        theme(axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.text.x  = element_text(size=18),
              axis.text.y = element_text(size = 18),
              plot.title = element_text(size=22),
              strip.text.x =  element_text(size=18))
              
              
    } else if (type == "h") {
      
      ggplot(Qprof, aes(x = pos.cM, xend = pos.cM, y = 0, yend = log10pval,
                        group = chr)) + geom_segment() +
        facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
        geom_hline(yintercept = threshold, colour = "red") + 
        theme_bw() + xlab("position [cM]") + ylab("-log10(p.val)") + 
        ggtitle(main) +
        theme(axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.text.x  = element_text(size=18),
              axis.text.y = element_text(size = 18),
              plot.title = element_text(size=22),
              strip.text.x =  element_text(size=18))
              
    }
    
  } else { # QTL info given
    
    if (type == "l") {
      
      ggplot(Qprof, aes(x = pos.cM, y = log10pval, group = chr)) + geom_line() + 
        facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
        geom_vline(aes(xintercept = pos.cM), pos.Q, linetype = "longdash",
                   colour = "black") +
        geom_hline(yintercept = threshold, colour = "red") + theme_bw() +
        xlab("position [cM]") + ylab("-log10(p.val)") + 
        ggtitle(main) + 
        theme(axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.text.x  = element_text(size=18),
              axis.text.y = element_text(size = 18),
              plot.title = element_text(size=22),
              strip.text.x =  element_text(size=18))
      
      
    } else if (type == "h") {
      
      ggplot(Qprof, aes(x = pos.cM, xend = pos.cM, y = 0, yend = log10pval,
                        group = chr)) + geom_segment() +
        facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
        geom_vline(aes(xintercept = pos.cM), pos.Q, linetype = "longdash",
                   colour = "black") +
        geom_hline(yintercept = threshold, colour = "red") + 
        theme_bw() + xlab("position [cM]") + ylab("-log10(p.val)") + 
        ggtitle(main) +
        theme(axis.title.x = element_text(size=18),
              axis.title.y = element_text(size=18),
              axis.text.x  = element_text(size=18),
              axis.text.y = element_text(size = 18),
              plot.title = element_text(size=22),
              strip.text.x =  element_text(size=18))
              
    }
    
  }
  
}