###################
# QTL_effect_plot #
###################

#'
#' Plot QTL effect accross model
#' 
#' Function to visualise and compare the effects of QTLs accross the different
#' possible models for the QTL effect (cross-specific, parental, ancestral,
#' bi-allelic).
#' 
#' The function calculates the effects of the different QTL positions in the
#' four possible models (cross-specific, parental, ancsetral, bi-allelic). Then
#' the effects are ploted using buble. The size of the buble is proportional to
#' the significance of the effect. The intensity of the colour is proportional
#' to the size of the effect. The red colour means that the allele has a negative
#' effect on the trait and the blue that it has a positive effect.
#' 
#' The effects are presented using the crosses as reference. For each cross,
#' the parent A and B are given. The effect of the cross-specific model are
#' represented per parent. The two parents have the same contribution because
#' the effect are estimated within cross contrasting the two parent effects.
#' One parent has a positive contribution, the other the same negative
#' contribution.
#' 
#' Then the effect of the parental model are plotted. By default (sum_zero = TRUE)
#' these effects are estimated using a sum to zero constraint. This means that within a
#' connected part the sum of the parental allelic effect is zero. These effects
#' are plotted using the parent per cross as reference. A colour index indicate
#' which parent are the shared between crosses.
#' 
#' The ancestral effect are also estimated using a sum to zero constraint by
#' default. The ancestral effect are translated at the parental level. Parents
#' belonging to the same ancestral group get the same effect.
#' 
#' Finally, the bi-allelic effect are also represented. The bi-allelic effects
#' are also estimated using a sum to zero constraint.
#' 
#' An alternative is to estimate the QTL effect of the parental, ancestral,
#' and bi-allelic model using a reference allele set to 0 (sum_zero = FALSE).
#' In that case, the most frequent parental, ancestral and SNP allele will
#' be used as reference. For the parental and the ancestral model it is also
#' possible to fix the reference allele using argument ref.par.
#' 
#' @param mppData An object of class \code{mppData}.
#' 
#' @param trait \code{Numerical} or \code{character} indicator to specify which
#' trait of the \code{mppData} object should be used. Default = 1.
#' 
#' @param QTL Object of class \code{QTLlist} representing a list of
#' selected position obtained with the function \code{\link{QTL_select}} or
#' vector of \code{character} marker positions names.
#' 
#' @param sum_zero Optional \code{Logical} value specifying if the QTL effect of
#' a parental or an ancestral model should be caculated using the sum to zero
#' constraint. Default = TRUE.
#' 
#' @param ref.par Optional \code{Character} expression defining the parental
#' allele that will be used as reference for the parental model. For the
#' ancestral model, the ancestral class containing the reference parent will be
#' set as reference. \strong{This option can only be used if the MPP design is
#' composed of a unique connected part}. Default = NULL.
#' 
#' @param trait.lab \code{Character} expression indicating the name of the trait.
#' Default = 'trait'.
#' 
#' @param text.size \code{Numeric} value specifying the size of graph axis text
#' elements. Default = 16.
#' 
#' @param output.loc Path to a folder where the results will be saved.
#' By default the function uses the current working directory.
#' 
#' @param plot_label \code{Character} for the mane of the plot. plot_label.pdf
#' will be saved at output.loc. Default = 'Qeff_tab'.
#' 
#' @return Return:
#' 
#' \item{Qeff_tab}{\code{List} with one element per QTL with the QTL effects
#' and p-values for the different models that where plotted.}
#' 
#' A file plot_label.pdf saved at output.loc containing the plot of the QTL
#' effects.
#' 
#' @author Vincent Garin
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(mppData)
#' 
#' SIM <- mpp_SIM(mppData = mppData)
#' QTL <- QTL_select(SIM)
#' 
#' # Specify a location where your results will be saved
#' my.loc <- "C:/.../..."
#' 
#' MQE <- QTL_effect_plot(mppData = mppData, QTL = QTL, output.loc = my.loc)
#'                  
#' }
#' 
#' @export

QTL_effect_plot <- function(mppData, trait = 1, QTL, sum_zero = TRUE,
                            ref.par = NULL, trait.lab = 'trait', text.size = 16,
                            output.loc = getwd(), plot_label = 'Qeff_tab'){
  
  if(!is.null(ref.par)){ sum_zero = FALSE }
  
  if(is.character(QTL)){
    
    Q_info <- mppData$map[mppData$map[, 1] %in% QTL,]
    
  } else {Q_info <- QTL}
  
  n_QTL <- dim(Q_info)[1]
  
  Q_title <- vector(mode = "list", length = n_QTL)
  
  for(i in 1:n_QTL){
    
    Qi <- paste0("Q", i)
    mki <- paste("mk:", Q_info[i, 1])
    chri <- paste("chr =", Q_info[i, 2])
    posi <- paste("pos =", round(Q_info[i, 4], 2), "cM")
    
    Q_title[[i]] <- paste(trait.lab, Qi, mki, chri, posi, sep = "; ")
    
  }
  
  # Estimation of the QTL effects of the given list for the 4 different models
  
  QEff_cr <- QTL_gen_effects(mppData = mppData, trait = trait, QTL = QTL,
                            Q.eff = 'cr')
  
  QEff_par <- QTL_gen_effects(mppData = mppData, trait = trait, QTL = QTL,
                             Q.eff = 'par', ref.par = ref.par,
                             sum_zero = sum_zero)
  
  QEff_anc <- QTL_gen_effects(mppData = mppData, trait = trait, QTL = QTL,
                             Q.eff = 'anc', ref.par = ref.par,
                             sum_zero = sum_zero)
  
  if(sum_zero){
    
    QEff_biall <- QTL_gen_effects(mppData = mppData, trait = trait, QTL = QTL,
                                 Q.eff = 'biall')
    
    # recalculate to get only the single beta values for the minor allele
    
    mppData2 <- mppData
    mppData2$geno.par <- NULL
    mppData2$geno.par.clu <- NULL
    
    QEff_biall2 <- QTL_gen_effects(mppData = mppData2, trait = trait, QTL = QTL,
                                  Q.eff = 'biall')
    
    rm(mppData2)
    
    # Need to convert the obtained result into coefficient obtained with sum to
    # zero constraint.
    
    Q_pos <- which(mppData$geno.par[, 1] %in% Q_info[, 1])
    
    for(i in 1:n_QTL){
      
      all_i <- mppData$allele.ref[, Q_pos[i], drop = FALSE]
      B_i <- QEff_biall2$Qeff[[i]]$Effect
      pval_i <- QEff_biall2$Qeff[[i]]$`p-value`
      all_sc <- c(B_i, -B_i, 0, 0)
      pval_sc <- c(pval_i, pval_i, 1, 1)
      names(all_sc) <- names(pval_sc) <- all_i
      Eff_i <- all_sc[QEff_biall$Qeff[[i]]$Par.all]
      QEff_biall$Qeff[[i]]$Effect <- Eff_i
      pval_i <- pval_sc[QEff_biall$Qeff[[i]]$Par.all]
      QEff_biall$Qeff[[i]]$`p-value` <- pval_i
      
    }
    
    
  } else {
    
    QEff_biall <- QTL_gen_effects(mppData = mppData, trait = trait, QTL = QTL,
                                 Q.eff = 'biall')
    
  }
  
  
  
  
  # Form the table of QTL effects
  ###############################
  
  n_cr <- mppData$n.cr
  
  cr_ind <- mppData$par.per.cross[, 1]
  cr_ind <- rep(cr_ind, each = 2)
  
  par_A <- mppData$par.per.cross[, 2]
  par_B <- mppData$par.per.cross[, 3]
  
  par_ind <- rep('', 2*n_cr)
  par_ind[seq(from = 1, to = (2*n_cr), by = 2)] <- par_A
  par_ind[seq(from = 2, to = (2*n_cr), by = 2)] <- par_B
  
  Qeff_tab <- vector(mode = 'list', length = n_QTL)
  
  for(i in 1:n_QTL){
    
    # cr
    
    add_ind <- (QEff_cr$Qeff[[i]]$Add.parent == mppData$par.per.cross[, 3]) * 1
    
    add_ind2 <- c()
    
    for(j in 1:n_cr){
      
      if(add_ind[j] == 0 || is.na(add_ind[j])){
        
        add_ind2 <- c(add_ind2, c(1, -1))
        
      } else {add_ind2 <- c(add_ind2, c(-1, 1)) }
      
    }
    
    add_eff_cr <- rep(QEff_cr$Qeff[[i]]$Effect, each = 2)
    add_eff_cr <- add_eff_cr * add_ind2
    p_val_cr <- rep(QEff_cr$Qeff[[i]]$`p-value`, each = 2)
    
    # par
    
    add_eff_par <- QEff_par$Qeff[[i]]$Effect
    p_val_par <- QEff_par$Qeff[[i]]$`p-value`
    names(add_eff_par) <- names(p_val_par) <- rownames(QEff_par$Qeff[[i]])
    
    add_eff_par <- add_eff_par[par_ind]
    p_val_par <- p_val_par[par_ind]
    
    # anc
    
    add_eff_anc <- QEff_anc$Qeff[[i]]$Effect
    p_val_anc <- QEff_anc$Qeff[[i]]$`p-value`
    names(add_eff_anc) <- names(p_val_anc) <- rownames(QEff_anc$Qeff[[i]])
    
    add_eff_anc <- add_eff_anc[par_ind]
    p_val_anc <- p_val_anc[par_ind]
    
    # biall
    
    add_eff_biall <- QEff_biall$Qeff[[i]]$Effect
    p_val_biall <- QEff_biall$Qeff[[i]]$`p-value`
    names(add_eff_biall) <- names(p_val_biall) <- rownames(QEff_biall$Qeff[[i]])
    
    add_eff_biall <- add_eff_biall[par_ind]
    p_val_biall <- p_val_biall[par_ind]
    
    
    df <- data.frame(cr_ind, par_ind, add_eff_cr, p_val_cr, add_eff_par, p_val_par,
                     add_eff_anc, p_val_anc, add_eff_biall, p_val_biall,
                     stringsAsFactors = FALSE)
    
    Qeff_tab[[i]] <- df
    
  }
  
  # form the colour legend
  ########################
  
  # form colour codes for parents ancestor and bi-allelic group
  
  par_ref <- 1:mppData$n.par
  names(par_ref) <- mppData$parents
  
  par_code <- par_ref[par_ind]
  
  # ancestral code
  
  par_clu_QTL <- mppData$par.clu[Q_info[, 1], , drop = FALSE]
  
  anc_code <- c()
  
  for(i in 1:n_QTL){
    
    par_clu_QTLi <- par_clu_QTL[i, ]
    anc_code_i <- par_clu_QTLi[par_ind]
    
    anc_ref <- 1:length(unique(anc_code_i))
    names(anc_ref) <- as.character(unique(anc_code_i))
    
    anc_code_i <- anc_ref[as.character(anc_code_i)]
    
    anc_code <- cbind(anc_code, anc_code_i)
    
  }
  
  # Bi-allelic code
  
  # bi_code <- c()
  # 
  # gen_par <- mppData$geno.par
  # 
  # SNP_QTL <- gen_par[gen_par[, 1] %in% Q_info[, 1], 5:dim(gen_par)[2], ]
  # 
  # SNP_QTL_ref <- mppData$allele.ref[, Q_info[, 1]]
  # 
  # 
  # for(i in 1:n_QTL){
  #   
  #   SNP_QTLi <- unlist(SNP_QTL[i, ])
  #   
  #   bi_code_i <- SNP_QTLi[par_ind]
  #   
  #   SNP_QTL_refi <- SNP_QTL_ref[, i]
  #   
  #   all_ref <- c(1, 3, 2, 2)
  #   names(all_ref) <- SNP_QTL_refi
  #   
  #   bi_code_i <- all_ref[bi_code_i]
  #   
  #   bi_code <- cbind(bi_code, bi_code_i)
  #   
  # }
  #
  # rownames(anc_code) <- rownames(bi_code) <-  par_ind
  # colnames(anc_code) <- colnames(bi_code) <- paste0('Q', 1:n_QTL)
  
  rownames(anc_code) <-  par_ind
  colnames(anc_code) <- paste0('Q', 1:n_QTL)
  
  # colour indicator plot
  #######################
  
  # Parent
  
  cr_par_ind <- paste("C" ,cr_ind, "P", par_ind, sep = "_")
  cr_par_ind <- rep(cr_par_ind, 4)
  y <- rep((2*n_cr):1, 4)
  y <- factor(y)
  x <- rep(1, length(cr_par_ind))
  model <- rep(c("cr", "par", "anc", "biall"), each = 2*n_cr)
  model <- factor(model, levels = c("cr", "par", "anc", "biall"))
  
  parAB_ind <- rep(c("; Pa: ", "; Pb: "), n_cr)
  
  y.names <- rev(paste0("Cr: " ,cr_ind, parAB_ind, par_ind))
  y.lab = ""
  
  y <- length(par_ind):1
  y <- factor(y)
  x <- rep(1, length(par_ind))
  z <- par_code
  z <- factor(z, levels = sort(unique(z)))
  ind <- rep("Par id", length(par_ind))
  n_col_par <- mppData$n.par
  
  data <- data.frame(x = x, y = y, z = z, ind = ind)
  
  # plot y axis and parent colour indicator
  
  col_par <- ggplot(data, aes(x, y, z = z)) +
    facet_wrap(nrow = 1, ~ ind, scales = "free_x") + geom_tile(aes(fill = z)) +
    scale_fill_manual(values = colorRampPalette( RColorBrewer::brewer.pal(9, "Set1"))(n_col_par)) +
    theme_bw() +  xlab("") + ylab(y.lab) + ggtitle("") +
    scale_y_discrete(labels = y.names)  + theme(legend.position='none') +
    theme(axis.text.y = element_text(size = text.size),
          axis.title.y = element_text(size = text.size),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = text.size),
          strip.text.x =  element_text(size=text.size),
          legend.title = element_text(size=(text.size)))
  
  # plot only the parent colour indicator
  
  # col_par <- ggplot(data, aes(x, y, z = z)) +
  #   facet_wrap(nrow = 1, ~ ind, scales = "free_x") + geom_tile(aes(fill = z)) +
  #   scale_fill_manual(values = colorRampPalette( RColorBrewer::brewer.pal(9, "Set1"))(colourCount)) +
  # theme(legend.position='none') +
  #   theme_bw() +  xlab("") + ylab(y.lab) + theme(legend.position='none') +
  #   theme(axis.text.y = element_blank(),
  #         axis.text.x = element_blank(),
  #         axis.title.x = element_text(size = text.size),
  #         strip.text.x =  element_text(size=text.size))
  
  # plot only legend
  
  # y_leg <- ggplot(data, aes(x, y, z = z)) + geom_blank() +
  #   ggtitle("") + scale_y_discrete(labels = y.names) + ylab("") +
  #   theme(axis.text.y = element_text(size = text.size),
  #         axis.line.x=element_blank(),
  #         axis.text.x=element_blank(),
  #         axis.ticks.x=element_blank(),
  #         axis.title.x=element_blank(),
  #         panel.grid.minor.x=element_blank(),
  #         panel.grid.major.x=element_blank())
  # 
  # plot_grid(y_leg, col_par, ncol = 2, rel_widths = c(0.5, 0.5))
  
  # Plot everything together
  ##########################
  
  # A loop will iterage over the ancestral and bi-allelic colour code
  # and put toghether these elements with the QTL effects
  
  pdf(file.path(output.loc, paste0(plot_label,".pdf")), height = 10, width = 16)
  
  for(i in 1:n_QTL){
    
    # ancestral colour
    
    z <- anc_code[, i]
    n_col_anc <- length(unique(z))
    z <- factor(z, levels = sort(unique(z)))
    ind <- rep("Anc id", length(par_ind))
    
    data <- data.frame(x = x, y = y, z = z, ind = ind)
    
    col_anc <- ggplot(data, aes(x, y, z = z)) +
      facet_wrap(nrow = 1, ~ ind, scales = "free_x") + geom_tile(aes(fill = z)) +
      scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_col_anc)) +
      theme(legend.position='none') +
      theme_bw() +  xlab("") + ylab(y.lab) + ggtitle("") +
      theme(legend.position='none') +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_text(size = text.size),
            strip.text.x =  element_text(size=text.size),
            legend.title = element_text(size=(text.size)))
    
    # bi-allelic colour
    
    # z <- bi_code[, i]
    # z <- factor(z, levels = sort(unique(z)))
    # cols <- c("1" = "black", "2" = "grey", "3" = "white")
    # ind <- rep("SNP id", length(par_ind))
    # 
    # data <- data.frame(x = x, y = y, z = z, ind = ind)
    # 
    # col_bi <- ggplot(data, aes(x, y, z = z)) +
    #   facet_wrap(nrow = 1, ~ ind, scales = "free_x") + geom_tile(aes(fill = z)) +
    #   scale_fill_manual(values = cols) + theme(legend.position='none') +
    #   theme_bw() +  xlab("") + ylab(y.lab) + theme(legend.position='none') +
    #   theme(axis.text.y = element_blank(),
    #         axis.text.x = element_blank(),
    #         axis.title.x = element_text(size = text.size),
    #         strip.text.x =  element_text(size=text.size))
    
    # QTL effect table
    
    graph_title <- Q_title[[i]]
    
    Q_eff_k <- Qeff_tab[[i]]
    
    effect_add <- c(Q_eff_k$add_eff_cr, Q_eff_k$add_eff_par,
                    Q_eff_k$add_eff_anc, Q_eff_k$add_eff_biall)
    pval <- c(Q_eff_k$p_val_cr, Q_eff_k$p_val_par, Q_eff_k$p_val_anc,
              Q_eff_k$p_val_biall)
    LOD <- -log10(pval)
    
    data <- data.frame(cr_par_ind, y, x, model, effect_add, pval, LOD,
                       stringsAsFactors = FALSE)
    
    p <- ggplot(data,aes(x=x,y=y, xmin=0.75, xmax=1.25 ,ymax=(2*n_cr),ymin=1,
                         size=LOD, color=effect_add)) +
      facet_wrap(nrow = 1, ~ model, scales = "free_x") +
      geom_point() + scale_size_continuous(range = c(0, 10)) +
      scale_colour_gradient2(low="red", high="blue", mid="white") +
      theme_bw() +  xlab("") + ylab("") + ggtitle(graph_title) +
      theme(axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            strip.text.x =  element_text(size=text.size),
            legend.title = element_text(size=(text.size)),
            legend.text = element_text(size=(text.size)))
    
    # print(plot_grid(col_par, col_anc, col_bi, p, ncol = 4,
    #                 rel_widths = c(0.2, 0.07, 0.07, 0.66)))
    
    print(cowplot::plot_grid(col_par, col_anc, p, ncol = 3,
                    rel_widths = c(0.23, 0.07, 0.7)))
    
  }
  
  dev.off()
  
  return(Qeff_tab)
  
}