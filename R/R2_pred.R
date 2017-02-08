###########
# R2_pred #
###########

# function to compute the predicted R squared in a validation set

# y.vs: observed trait value in the validation set

# B.ts: list of QTL genetic effects estimated in the training set

# Q.list: list of QTL incidence matrices

# within.cross: indicate if the prediction must be computed within cross


R2_pred <- function(mppData.vs, B.ts, Q.list, within.cross) {
  
  y.vs <- mppData.vs$trait
  
  proj.fct <- function(X, B){
    
    ############# option here to detect outliers in Beta
    
    X[is.na(X)] <- 0; B[is.na(B)] <- 0 # replace missing values by 0
    X %*% B
    
  }
  
  y.pred <- mapply(proj.fct, X = Q.list, B = B.ts)
  y.pred <- rowSums(y.pred)
  
  dataset <- cbind(y.vs, y.pred)
  cross.ind <- mppData.vs$cross.ind
  
  index <- complete.cases(dataset)
  dataset <- dataset[index, ]
  cross.ind <- cross.ind[index]
  
  if(within.cross){
    
    with.cross.cor <- function(x){
      if((length(unique(x[, 1])) == 1) || (length(unique(x[, 2])) == 1)){ NA
      } else {cor(x[, 1], x[, 2])^2}
    }
    
    dataset.cr <- split(x = dataset,
                        f = factor(cross.ind, levels = unique(cross.ind)))
    R2 <- lapply(X = dataset.cr, FUN =  with.cross.cor)
    R2 <- 100 * (mean(unlist(R2), na.rm = TRUE))
    
  } else { R2 <-  100 * (cor(dataset[, 1], dataset[, 2])^2) }
  
  return(R2)
  
}