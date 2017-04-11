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
#' effects. White colour means no significant effect.
#' 
#' For a cross-specific QTL profile (\code{Q.eff = "cr"}): Red color means
#' that the allele coming from parent A(1) increases the phenotypic value and
#' parent B(2) decreases it and blue that parent A(1) decreases the trait and
#' parent B(2) increases it.
#' 
#' For a parental model (\code{Q.eff = "par"}): Effect are defined within
#' connected part (indicated as ci). Red (Blue) colour means a signicative
#' negative (positive) effect with respect to the reference of the connected
#' part. The reference alleles are defined as the most used one regarding first
#' the number of cross where it segregates and second the number of genotypes
#' where it potentially segregates.
#' 
#' For an ancestral model the the parent with a red (blue) colour receive an
#' allele that has a positive (negative) effect with respect to the reference
#' ancestral allele of the connected part. Since the results of the ancetral
#' clustering potentially change at each position. It is not possible to
#' indicate which allele are the reference alleles. But the reference allele
#' within each connected part are defined as for the parental model as the one
#' that segregate the most. Therefore interpretation
#' of the genetic effect plot should be done with caution. In that case, the
#' plot should be taken as a rough indication of the signal distribution.
#' 
#' The colour intensity increase with the significance of the effect.
#' There are five colour intensities according to the p-value of the QTL
#' effect: 0.05<p-val<0.01; 0.01<p-val<0.001; 0.001<p-val<0.0001;
#' 0.0001<p-val<0.00001 and p-val< 0.00001.
#' 
#' The dashed lines indicate the position that have been introduce in the
#' argument \code{QTL}.
#' 
#' @param mppData An object of class \code{mppData}.
#' See \code{\link{mppData_form}} for details.
#' 
#' @param Qprof Object of class \code{QTLprof} returned by the function
#' \code{\link{mpp_SIM}} or \code{\link{mpp_CIM}} with argument
#' \code{est.gen.eff = TRUE}.
#' 
#' @param Q.eff \code{Character} expression indicating the assumption concerning
#' the QTL effect: 1) "cr" for cross-specific effects; 2) "par" parental
#' effects; 3) "anc" for an ancestral effects.
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
#' # without QTL positions
#' plot_genEffects(mppData = USNAM_mppData, Qprof = SIM, Q.eff = "cr") 
#' 
#' # with QTL positions
#' plot_genEffects(mppData = USNAM_mppData, Qprof = SIM, Q.eff = "cr",
#' QTL = QTL) 
#' 
#' @export
#' 


plot_genEffects <- function(mppData, Qprof, Q.eff, QTL = NULL,
                            main = "QTL genetic effects plot")
{
  
  # 1. check data format
  ######################
  
  stopifnot(inherits(Qprof, "QTLprof"))
  
  if(!(Q.eff %in% c("cr", "par", "anc"))){
    
    stop("The Q.eff argument must take value: 'cr', 'par', 'anc'.")
    
  }
  
  
  n.eff <- dim(Qprof)[2] - 5
  
  if(n.eff == 0) {
    
    stop("The Qprof object does not contain any QTL p-value information.
         It was probably not obtained using est.gen.eff = TRUE")
    
  }
  
  # 2. order columns within connected parts
  #########################################
  
  if((Q.eff == "par") || (Q.eff == "anc")){
    
    # determine the connected parts
    
    con.part <- design_connectedness(par.per.cross = mppData$par.per.cross,
                                     plot.des = FALSE)
    
    allele_order <- c()
    pval <- data.frame(row.names = 1:dim(Qprof)[1])
    
    for(i in seq_along(con.part)){
      
      con.part_i <- con.part[[i]]
      
      # subset results of the connected part
      
      pval_i <- Qprof[, con.part_i]
      
      # order according to number of reference values
      
      all.ref <- apply(X = pval_i, MARGIN = 2,
                       FUN = function(x) sum(x == 1))
      
      pval <- cbind.data.frame(pval, pval_i[, names(sort(all.ref))])
      
      allele_ord_i <- names(sort(all.ref))
      
      allele_ord_i <- c(paste(allele_ord_i, paste0("(c", i,")"), sep = "\n"))
      
      allele_order <- c(allele_order, allele_ord_i)
      
    }
    
    # Rename Qprof
    
    Qprof <- cbind(Qprof[, 1:5], pval)
    
    
    y.names <- allele_order
    
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
  
  chr <- rep(Qprof$chr, n.eff)
  
  ### 2.4 genetic map positions in cM with width between two positions.
  
  x <- rep(Qprof$pos.cM, n.eff)
  
  w <- tapply(X = Qprof$pos.cM, INDEX = Qprof$chr,
              FUN = function(x) c(diff(x), 1))
  w <- unlist(w)
  w <- rep(w, n.eff)
  
  pos.cM <- Qprof$pos.cM
  
  y_lab <- "parents"
  
  ### 2.5 legend for the y-axis (cross or parents)
  
  if(Q.eff == "cr") {
    
    cross.names <- unique(mppData$cross.ind)
    par.cross.names <- paste0("(", mppData$par.per.cross[, 2], 
                              "x", mppData$par.per.cross[, 3], ")")
    
    y.names <- paste(cross.names, par.cross.names, sep = "\n")
    
    y_lab <- "crosses"
    
  } 
  
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
      theme_bw() + xlab("position [cM]") + ylab(y_lab) + 
      scale_y_discrete(labels = y.names) + ggtitle(main) +
      theme(axis.title.x = element_text(size=18),
            axis.title.y = element_text(size=18),
            axis.text.x  = element_text(size=18),
            axis.text.y = element_text(size = 18),
            plot.title = element_text(size=22),
            strip.text.x =  element_text(size=18),
            legend.title = element_text(size=16),
            legend.text = element_text(size=16))
    
  } else { # QTL position given
    
    pl <- ggplot(data, aes(x, y, z = z))
    pl + geom_tile(aes(fill = z, width = w)) +
      facet_wrap(nrow = 1, ~ chr, scales = "free_x") +
      scale_fill_gradient2(limits = c(-5, 5), low = "red", mid = "white",
                           high = "blue") +
      geom_vline(aes(xintercept = pos.cM), pos.Q, linetype = "longdash") +
      theme_bw() + xlab("position [cM]") + ylab(y_lab) +
      scale_y_discrete(labels = y.names) + ggtitle(main) +
      theme(axis.title.x = element_text(size=18),
            axis.title.y = element_text(size=18),
            axis.text.x  = element_text(size=18),
            axis.text.y = element_text(size = 18),
            plot.title = element_text(size=22),
            strip.text.x =  element_text(size=18),
            legend.title = element_text(size=16),
            legend.text = element_text(size=16))
    
  }
  
}