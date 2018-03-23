#############################################################################
###             Get the probability of each leafs of the tree           #####
### Based on the tree structure, the probability can calculate on each  #####
### node by the multiplicate along the tree                             #####
### prob is the probability on each branch                              #####
#############################################################################

#-----------------------------------------------------------------------#####

Prob.mult = function(p, tree, prob){
  
  prob_mult <- rep(0, p)
  
  for(t in 1:(p + tree$Nnode - 1)){
    index <- which(tree$edge[t, 1] == tree$edge[, 2])
    
    if(length(index) == 0){
      prob_mult[t] <- prob[t]
    }
    else{
      prob_mult[t] <- prob[t] * prob_mult[index]
    }
  }
  prob_leaf <- prob_mult[tree$edge[, 2] <= p]
  
  if(abs(sum(prob_leaf) - 1) < 1e-20){
    result <- data.frame(prob_leaf)
    return(result)
  }
  else{
    stop("Error: Sum of leafs is not equal 1.")
  }
}

##############################################################################
