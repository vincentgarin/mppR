#################################
# function for clusthaplo plots #
#################################

.use.deprecated.draw.chromosomes <- F

plot.clusthaplo.haplotypes <- function(x, ...) {
  sup <- list(...)
  if(is.null(sup$ord)) {
    sup$ord <- 1:ncol(x)
  }
  if(.use.deprecated.draw.chromosomes) {
    if(is.null(sup$names)) {
      sup$names <- attr(x, 'ind.names')[sup$ord]
    }
    if(is.null(sup$loci)) {
      sup$loci <- attr(x, 'loci')
    }
    sup$cl <- .build.cliq.struc(.cliques(x))
    do.call(.draw.chromosome.OBSOLETE, sup)
    title(xlab="Map (cM)")
  } else {
    sup$tc <- x
    sup$step.size <- max(diff(attr(x, 'loci')))
    do.call(.draw.chrom, sup)
  }
}

.DSAT <- function(cl, R=cl$repr) {
  graph <- .cliq.graph(cl, R)
  color.vec <- rep(NA, ncol(graph))
  singlz <- .singletons(cl)
  color.vec[singlz] <- 0  # forcefully set this color.
  
  degrees <- colSums(graph)
  
  dsat.value <- function(v) {
    color.nei <- unique(na.omit(color.vec[graph[, v]]))
    color.nei <- color.nei[color.nei != 0]
    if (length(color.nei) == 0) {
      return(degrees[v])
    } else {
      return(sum(color.nei != 0))  # modified saturation
    }
  }
  
  pick.vertex <- function() {
    candidates <- is.na(color.vec)
    if(!any(candidates)) {
      return(NULL)
    }
    candidates <- which(candidates)
    dsat.vec <- sapply(candidates, dsat.value)
    strongest <- candidates[dsat.vec == max(dsat.vec)]
    #cat("candidate vertices", candidates, ".DSAT", max(dsat.vec), "\n")
    if (length(strongest) > 1) {
      ret <- strongest[which.max(degrees[strongest])]
    } else {
      ret <- strongest
    }
    #cat("selecting vertex", ret, "\n")
    ret
  }
  
  pick.color <- function(v) {
    color.nei <- unique(na.omit(color.vec[graph[, v]]))
    if(max(c(0, color.nei)) == 0) {
      return(1)
    }
    candidates <- 1:(1+length(color.nei))
    ret <- which(!(candidates %in% color.nei))[1]
    #cat("selecting color", ret, "\n")
    ret
  }
  
  while (!is.null(v <- pick.vertex())) {
    color.vec[v] <- pick.color(v)
    #cat("color[", v, "] = ", color.vec[v], "\n", sep="")
    #cat(color.vec, "\n")
  }
  color.vec
}

.draw.chrom <- function(tc, ord=1:ncol(tc), expand.cliques=F, step.size=1, from=min(attr(tc, 'loci')), to=max(attr(tc, 'loci')), background='gray40') {
  x <- NULL  # make Rcheck happy
  y <- NULL  # make Rcheck happy
  cl <- .build.cliq.struc(.cliques(tc))
  if(length(ord) != ncol(cl$repr)) {
    # need to rebuild cliques
    reord <- 1:ncol(cl$repr)
    reord[ord] <- 1:length(ord)
    sel <- reord[.select(.useable.index(cl), unique(c(ord, 1:ncol(cl$repr))))]
    repr <- matrix(sel[cl$repr], nrow=nrow(cl$repr), ncol=ncol(cl$repr))[, ord]
    cl <- .build.cliq.struc(.cliques(repr))
  }
  names <- attr(tc, "ind.names")
  if(expand.cliques == 'compare') {
    tmp <- .maximize.repr(cl)[, ord]
    R <- matrix(0, ncol=2 * ncol(tmp), nrow=nrow(tmp))
    R[, (1:ncol(tmp)) * 2] <- tmp
    R[, (1:ncol(tmp)) * 2 - 1] <- cl$repr
    tmp <- names
    #print(tmp)
    names[(1:length(tmp)) * 2] <- paste(tmp, ".m", sep="")
    names[(1:length(tmp)) * 2 - 1] <- tmp
    tmp <- ord
    ord[(1:length(tmp)) * 2] <- tmp * 2
    ord[(1:length(tmp)) * 2 - 1] <- tmp * 2 - 1
  } else if(expand.cliques) {
    R <- .maximize.repr(cl)[, ord]
  } else {
    R <- cl$repr[, ord]
  }
  
  .colors <- matrix(.DSAT(cl, R)[R], ncol=ncol(R))
  #print(dim(.colors))
  
  loci <- attr(tc, 'loci')
  
  to.keep <- loci <= to & loci >= from
  
  .colors <- .colors[to.keep, ]
  loci <- loci[to.keep]
  
  ydf <- data.frame(y=1:ncol(.colors))
  #geom.chromlines <- geom_hline(yintercept=1:ncol(.colors))
  
  cldf <- data.frame(x=rep(c(from, to), ncol(R)),
                     y=rep(1:ncol(R), each=2))
  geom.chromlines <- geom_line(aes(x=x, y=y, group=y),
                               data=cldf)
  breaks <- 1 + which(apply(diff(R)!=0, 1, any))
  
  xl <- xlim(from, to)
  yl <- #ylim(0, 1 + ncol(.colors))
    ylim(names)
  
  ncolors <- max(max(.colors) + 1, 2)
  
  .info.out("Using", ncolors, "colors.\n")
  
  ret <- (.gp(background=background)
          + xl + yl
          + xlab("Map (cM)")
          #+ ylab("Haplotypes")
          + ylab("")
          + geom.chromlines
  )
  if(length(breaks) > 0) {
    geom.loci <- geom_vline(xintercept=loci[-breaks],
                            colour='#404040C0', size=.1)
    geom.chrom <- .rect.mat(loci, .colors, step.size)
    geom.breaks <- geom_vline(xintercept=loci[breaks],
                              colour='#000000C0', size=.2)
    ret <- ret + geom.loci + geom.breaks + geom.chrom
  } else {
    geom.loci <- geom_vline(xintercept=loci,
                            colour='#404040C0', size=.1)
    ret <- ret + geom.loci
  }
  ret + scale_colour_manual(values=rainbow(ncolors))
}

.info.out <- function(...) .dumper(.verbose, ...)

.dumper <- function(flag, ...) {
  if (flag) {
    argv <- list(...)
    if (length(argv) == 1 && is.list(argv[[1]])) {
      print(argv[[1]])
    } else {
      cat(...)
    }
  }
}

.verbose <- F

.gp <- function(background='gray80') {
  (ggplot()
   + theme_bw()
   + theme(panel.background=element_rect(fill=background, size=0),
           panel.grid=element_blank(),
           panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),
           panel.spacing=element_blank(),
           axis.ticks.length=unit(.1, "cm")
   )
  )
}

.rect.mat <- function(locus.vec, col.mat, step=1) {
  x <- NULL  # make Rcheck happy
  y <- NULL  # make Rcheck happy
  height <- NULL  # make Rcheck happy
  width <- c(diff(locus.vec), step)
  make.height <- function(col.vec) {
    height <- rep(.8, length(width))
    height[col.vec==0] <- .01
    height
  }
  xvec <- rep(locus.vec + width * .5, ncol(col.mat))
  wvec <- rep(width, ncol(col.mat))
  hvec <- as.vector(apply(col.mat, 2, make.height))
  yvec <- rep(1:ncol(col.mat), each=length(locus.vec))
  fillvec <- 1 + as.integer(as.vector(col.mat))
  d <- data.frame(x=xvec, y=yvec, width=wvec, height=hvec)
  geom_tile(aes(x=x, y=y, width=width, height=height), fill=fillvec, data=d)
}

.build.cliq.struc <- function(cliqlist) {
  ci <- .clique.index(cliqlist)
  cliqlist <- lapply(cliqlist, function(c) { sapply(c, .index.of.clique, ci) })
  list(index=ci,
       cliques=cliqlist,
       n.ind=max(sapply(ci, max)),
       repr=.cliq.table(ci, cliqlist))
}

.clique.index <- function(cliqlist) {
  unique(unlist(cliqlist, recursive=F))
}

.cliques <- function(effet) {
  one_clique <- function(l, i) which(effet[l, ] == i)
  clique_list <- function(l) {
    x <- list()
    xi <- 1
    #for(i in 1:dim(effet)[2]) {
    for(i in sort(unique(effet[l, ]))) {
      o <- one_clique(l, i)
      if(length(o) > 0) {
        #cat(l, i, ":", o, "\n")
        x[[xi]] <- o
        xi <- xi + 1
      }
    }
    x
  }  
  lapply(1:dim(effet)[1], clique_list)
}

.index.of.clique <- function(cliq, cliqindex) {
  for(i in 1:length(cliqindex)) {
    ci <- cliqindex[[i]]
    if(length(ci) == length(cliq) && all(cliq == ci)) {
      return(i)
    }
  }
}

.cliq.table <- function(idx, cli) {
  n.ind <- max(sapply(idx, max))
  ret <- matrix(ncol=n.ind, nrow=length(cli))
  for(l in 1:length(cli)) {
    for(c in cli[[l]]) {
      ret[l, idx[[c]]] <- c
    }
  }
  ret
}

.cliq.graph <- function(cl, R=cl$repr) {
  graph <- matrix(F, ncol=length(cl$index), nrow=length(cl$index))
  #R <- cl$repr
  
  cliqz <- unique(R[1, ])
  
  graph[cliqz, cliqz] <- T
  
  for(pos in 2:nrow(R)) {
    new.cliqz <- unique(R[pos, ])
    graph[new.cliqz, new.cliqz] <- T
    graph[new.cliqz, cliqz] <- T
    graph[cliqz, new.cliqz] <- T
    cliqz <- new.cliqz
  }
  graph[col(graph) == row(graph)] <- F
  graph
}

.singletons <- function(cl) {
  unlist(lapply(1:length(cl$index),
                function(cli.num) {
                  if (length(cl$index[[cli.num]]) == 1) {
                    cli.num
                  }
                }))
}