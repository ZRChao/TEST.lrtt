#############################################################################
###            Tree.Ratio based on the tree structure                   #####
### simulate the tree structure based on the pacakge ape                #####
### which need define the OTU numbers p and tip.labed is 1:p            #####
### all the edge.length equals to 1                                     #####
#############################################################################

#-----------------------------------------------------------------------#####


Tree.ratio = function(p, tree, taxa.index, all.tab, group = c(1:(2*N))){
  result <- list()
  taxa_pv <- rep(0, ncol(taxa.index))
  names(taxa_pv) < names(sort(colSums(taxa.index)))
  
  names(otu_dif) <- as.character(1:p)
  taxa_leafs <- unique(sort(colSums(taxa.index)))
  
  otu_dif <- rep(F, p)
  taxa.dif <- c()
  n = 0
  for(t in taxa_leafs){
    taxa_names <- names(which(colSums(taxa.index) == t))
    for(j in 1:length(taxa_names)){
      n <- n + 1
      parent <- taxa_names[j]
      childs <- tree$edge[which(tree$edge[, 1] == parent), 2]

      child1 <- all.tab[, colnames(all.tab) == as.character(childs[1])]
      child2 <- all.tab[, colnames(all.tab) == as.character(childs[2])]

      ratio = log(child1 + 1e-20) - log(child2 + 1e-20)
      taxa_pv[n] <- pv <- ifelse(length(unique(ratio)) == 1, 1,
                                 t.test(ratio[group==1], ratio[group==2])$p.value)
                                 
      if(pv <= 0.05/length(taxa_names)){
        #temp = all_table[, colnames(all_table) == parent]
        all.tab[, colnames(all.tab) == parent] <- 1
        childs <- which(taxa.index[, colnames(taxa.index) == parent] == 1)
        taxa.dif <- append(taxa.dif, parent)
        
        if(t != p){
          otu_dif[childs] <- T
        }
        else{
          break
        }
      }
      else{
        all.tab[, which(colnames(all.tab) == parent)] <- child1 + child2
      }
    }
  }
  result$taxa.pvalue <- taxa_pv
  result$otu.dif <- otu_dif
  result$alltab <- all.tab
  result$taxa.dif <- taxa.dif
  return(result)
}

###########################################################################################

Tree.ratio.back = function(p, tree.ratio, taxa.index, otutab, group=c(1:(2*N))){
    
    difset <- which(tree_ratio$otu.dif == T)
    if(length(difset) == 0){
        difset.detected <- NULL
    }
    else{
        nondifset <- setdiff(1:p, difset)
        
        if(length(nondifset) == 0){
            ratio_sum <- log(otutab[, difset] + 1e-20) - log(rowSums(otutab) + 1e-20)
        }
        
        if(length(nondifset) == 1 ){
            ratio_sum <- log(otutab[, difset] + 1e-20) - log(otutab[, nondifset] + 1e-20)
        }
        
        else{
            ratio_sum <- log(otutab[, difset] + 1e-20) - log(rowSums(otutab[, nondifset]) + 1e-20)
        }
        difset_sumpv <-  apply(as.matrix(ratio_sum), 2, function(x)
            return(ifelse(length(unique(x)) == 1, 1, t.test(x[group==1], x[group==2])$p.value)))
            
        difset_sumpv.adj <- p.adjust(difset_sumpv, method = "fdr", n = length(difset_sumpv))
        difset.detected <- names(which(difset_sumpv.adj <= 0.05))
    }
    return(difset.detected)
}

#############################################################################################

