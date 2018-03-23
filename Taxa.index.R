#############################################################################
###            Get the relationship of internal nodes with leaf nodes    ####
### we can transfrom the tree structure defined by the edge pairs to     ####
### one matrix we call it taxa.index and each row is one leaf relation-  ####
### -ship with the internal nodes : 1 means it can connected 0 is not    ####
###    This function have the same function with caper::clade.members.list ##
###                                                                      ####
#############################################################################

#-----------------------------------------------------------------------#####

Taxa.index = function(p, phy_tree){
  taxa_index <- matrix(0, p, phy_tree$Nnode)
  colnames(taxa_index) <- as.character((p+1):(p+phy_tree$Nnode))
  for(i in 1:p){
   temp <- i
    for(j in 1:1000){
      index <- which(phy_tree$edge[, 2] == temp)
      if(length(index) != 0){
        parent <- phy_tree$edge[index, 1]
        taxa_index[i, which(colnames(taxa_index) == as.character(parent))] <- 1
        temp <- parent
      }
     else{
        break
      }
    }
  }
  return(taxa_index)
}

# so it is easy get the internal nodes counts
# taxa = otu_table%*%taxa_index
# colnames(taxa) <- as.character((p+1):(p+phy_tree$Nnode))
# all_table = cbind(otu_table, taxa)

#############################################################################
