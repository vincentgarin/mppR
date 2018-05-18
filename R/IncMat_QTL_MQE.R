##################
# IncMat_QTL_MQE #
##################

# function to produce different type of QTL incidence matrices

IncMat_QTL_MQE <- function(x, mppData, mppData_bi, Q.eff, par.clu,
                           cross.mat, par.mat, order.MAF){
  
  if(Q.eff == "biall") {
    
    inc_mat_QTL(x = x, mppData = mppData_bi, Q.eff = Q.eff)
    
  } else {
    
    inc_mat_QTL(x = x, mppData = mppData, Q.eff = Q.eff, order.MAF = order.MAF)
    
  }
  
}