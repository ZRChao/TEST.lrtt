#############################################################################
###            Multinomial tree distribution simulation                 #####
###   Given one tree structure and sample depth, on each internal nodes #####
### it will follow multinomial distribution with probability vectors,   #####
### here, depth is depth is unifrom p*10, p*1000 both case and control  #####
### N is sample size for case and control,  prob_control, prob_case is  #####
### probability on each branch.                                         #####
###    The results is the count data on leafs and internal nodes which  #####
### colname is correspond to the tree structure                         #####
#############################################################################

#-----------------------------------------------------------------------#####
BIT.Sim = function(p, seed= 2, N = 50, tree, prob_control, prob_case){
  count <- matrix(NA, 2*N, p + tree$Nnode - 1)
  
  for(j in (p + 1) : (p + tree$Nnode)){
    prob_control[which(tree$edge[, 1] == j)][1] -> prob1
    prob_case[which(tree$edge[, 1] == j)][1] -> prob2 
    
    for(i in 1:N){
      if(j == (p+1)){
        set.seed(i*j*seed)
        parent1 <- round(runif(1, p*10, p*1000))
        set.seed(i*j*(seed+1))
        parent2 <- round(runif(1, p*10, p*1000))
      }
      else{
        parent1 <- count[i, which(tree$edge[, 2] == j)]
        parent2 <- count[i + N, which(tree$edge[, 2] == j)]
      }
      
      set.seed(i*j*seed)
      child11 <- rbinom(1, parent1, prob1)
      child12 <- parent1 - child11
      count[i, which(tree$edge[, 1] == j)] <- c(child11, child12)
      
      set.seed(i*j*(seed+1))
      child21 <- rbinom(1, parent2, prob2)
      child22 <- parent2 - child21
      count[i+N, which(tree$edge[, 1] == j)] <- c(child21, child22)
    }
  }
  colnames(count) <- as.character(phy_tree$edge[, 2])
  return(count)
}


#############################################################################
