#############################################################################
###       logistical normal multinomial distribution simulation         #####
###   Given one tree structure and sample depth, on leafs it will follow#####
### logistical normal multinomial distribution. Correspond to the tree, #####
### choose the different OTU with different mu in case and control, and #####
### then take logit transform as the probability on the leafs count     #####
### here, depth is depth is unifrom p*10, p*2000 both case and control  #####
### N is sample size for case and control.                              #####
###    The results is the count data on leafs and combine internal node #####
### counts which is the sum of his childs, the colname is correspond to #####
### the tree structure.                                                 #####
#############################################################################

#-----------------------------------------------------------------------#####

LNM.Sim = function(p, seed, N, taxa_index = taxa_index, dif = diff_leaf){
  lnm_otutab <- matrix(NA, 2*N, p)
  set.seed(seed)
  control_mu <- runif(p-1, -0.00001, 0.00001)
  
  # need to import mvtnorm
  for(y in 1:N){
    set.seed(seed*y+5)
    control_y <- rmvnorm(1, mean = control_mu, sigma = diag(length(control_mu)))
    control_z <- c(exp(control_y)/(sum(exp(control_y))+1), 1/(sum(exp(control_y))+1))
    set.seed(y + 1)
    control_M <- round(runif(1, p*10, p*2000))
    lnm_otutab[y, ] <- round(control_z * control_M)
  }
  
  case_mu <- control_mu
  diff_num <- length(dif)
  set.seed(p)
  dif_mu <- c(runif(round(diff_num/2) ,-4, -1),runif(diff_num -round(diff_num/2) ,1, 4))
  case_mu[dif] <- dif_mu
  
  for(x in 1:N){
    set.seed(seed*x+3)
    case_y <- rmvnorm(1, mean = case_mu, sigma = diag(length(case_mu)))
    case_z <- c(exp(case_y) / (sum(exp(case_y))+1), 1/(sum(exp(case_y))+1))
    set.seed(seed*x + 1)
    case_M <- round(runif(1, p*10, p*2000))
    set.seed(seed*x + 2)
    lnm_otutab[x+N, ] <- rmultinom(1, case_M, case_z)
  }
  colnames(lnm_otutab) <- as.character(1:p)
  taxa <- lnm_otutab%*%taxa_index
  lnm_alltab <- cbind(taxa, lnm_otutab)
  return(lnm_alltab)
}

#############################################################################
