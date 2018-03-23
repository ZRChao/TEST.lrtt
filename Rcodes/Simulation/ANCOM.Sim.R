#############################################################################
###            ANCOM: Poission distribution simulation                  #####
###   Each OTU on the leaf is sampled from poission distribution.       #####
###   Correspond to the tree structure, the differential OTU on the     #####
###   tree leafs have different mu for poission distribution otherwise  #####
###   follow the same mu. The results is the count data on leafs and    #####
###   colname is correspond to the tree structure. All the parameters   #####
###   set accroding to ANCOM simulation.                                #####
#############################################################################

#-----------------------------------------------------------------------#####

ANCOM.Sim = function(p, seed = 1, N, dif = diff_leaf){
  control <- case <- matrix(NA, N, p)
  
  for(i in 1:N){
    a <- 200
    for(o in 1:p){
      A <- c()
      B <- c()
      for(n in 1:N){
       set.seed((o + n)*i*seed)
        mu <- rgamma(1, a, 1)
        set.seed((o*n)*i*seed)
        A[n] <- rpois(1, mu)
      
        if(o %in% dif){
          ul <- 200
          set.seed((o + n)*i*seed + 2)
          u <- runif(1, ul, 3*ul/2)
          set.seed((o + n)*i*seed + 3)
          B[n] <- rpois(1, mu + u)
        }
        else{
          set.seed((o + n)*i*seed + 4)
           B[n] <- rpois(1, mu)
        }
      }
      case[, o] <- B
      control[, o] <- A
    }
  }  
  ancom_otutab <- rbind(case, control)
  colnames(ancom_otutab) <- as.character(1:p)
  return(ancom_otutab)
}

#############################################################################


