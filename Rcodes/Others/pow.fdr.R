#############################################################################
###            Calculate the power and fdr by adjust pvalue             #####
###                                                                     #####
#############################################################################

#-----------------------------------------------------------------------#####

pow.fdr = function(p.adj, dif = diff_leaf, cutf = 0.05){
    
  reject = which(p.adj <= cutf)
  TP = length(intersect(reject, dif))
  FP = length(reject) - TP
  
  power = TP/length(dif)
  fdr = ifelse(length(reject) == 0, 0, FP/length(reject))
  return(c(power, fdr))
}

#############################################################################
