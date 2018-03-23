#############################################################################
###       Dirichlet Multinomial tree distribution simulation            #####
###   Given one tree structure and sample depth, on each internal nodes #####
### it will follow dirichlet multinomial distribution with probability  #####
### here, depth is depth is unifrom p*10, p*1000 both case and control  #####
### N is sample size for case and control,case_pi and control_pi is the #####
### probability on each branch for dirichlet distribution.              #####
### the dispersion theta is ranges from 0.001-0.1 on lower internal node#####
###    The results is the count data on leafs and internal nodes which  #####
### colname is correspond to the tree structure                         #####
#############################################################################

#-----------------------------------------------------------------------#####

DTM.Sim = function(p, seed = 1, N, tree, taxa_index = taxa_index, control_pi,
case_pi, theta = 0.1){
  
  dtm_alltab <- matrix(NA, 2*N, p + tree$Nnode - 1)
  colnames(dtm_alltab) <- as.character(tree$edge[, 2])
 
 for(j in (p + 1) : (p + tree$Nnode)){
    control_pi[which(tree$edge[, 1] == j)] -> pi0
    case_pi[which(tree$edge[, 1] == j)] -> pi1 
    theta0 <- ifelse(sum(taxa_index[, colnames(taxa_index) == as.character(j)]) > p/5,
        1e-8, theta)
        
    for(i in 1:N){
      if(j == (p+1)){
        set.seed(i*j*seed)
        parent1 <- round(runif(1, p*10, p*1000))
        set.seed(i*j*(seed+1))
        parent2 <- round(runif(1, p*10, p*1000))
      }
      else{
        parent1 <- dtm_alltab[i, which(tree$edge[, 2] == j)]
        parent2 <- dtm_alltab[i + N, which(tree$edge[, 2] == j)]
      }
      set.seed(i*j*seed)
      # need to import rdirmult packages
      dtm_alltab[i, which(tree$edge[, 1] == j)] <-
          simPop(J=1, K=2, n=parent1, pi = pi0, theta = theta0)$data
      set.seed(i*j*(seed+1))
      dtm_alltab[i+N, which(tree$edge[, 1] == j)] <-
          simPop(J=1, K=2, n=parent2, pi = pi1, theta = theta0)$data
    }
  }
  return(dtm_alltab)
}

#############################################################################
