#############################################################################
###            Simulate the probability on each branch                  #####
###                                                                     #####
#############################################################################

#-----------------------------------------------------------------------#####

# choose the internal node which will have different probability for his 
# children on different groups
set.seed(s)
diff_node = sample(c((p+1) : (p + phy_tree$Nnode))[which(colSums(taxa_index)<=p/10)], 5);diff_node


control_pi <- rep(0, phy_tree$Nnode)
for(j in (p + 1) : (p + phy_tree$Nnode)){
  set.seed(j)
  pi0 <- runif(1, 0.1, 0.9)
  control_pi[which(phy_tree$edge[, 1] == j)] <- c(pi0, 1 - pi0)
}

control_pi -> case_pi   ###case_probability
for(i in 1:length(diff_node)){
  parent <- diff_node[i]
  temp_diff <- which(phy_tree$edge[, 1] == parent)
  case_pi[temp_diff] <- control_pi[rev(temp_diff)]
}

tree_info <- cbind(phy_tree$edge, control_pi, case_pi)

## use the Prob.mult get the probability on leaf nodes and get the differential OTU

prob.leaf.g1 <- Prob.mult(p, phy_tree, case_pi)
prob.leaf.g2 <- Prob.mult(p, phy_tree, case_pi)
diff_leaf = which(prob.leaf.g1 != prob.leaf.g2)

phy_tree$diff_taxa <- diff_node
phy_tree$diff_otu <- diff_leaf
save(phy_tree, file = "phy.tree.rda")

# check is there any degeneration nodes

diff_otu = c()
for(i in 1:length(diff_node)){
    diff_child = which(taxa_index[, which(colnames(taxa_index) == as.character(diff_node[i]))] == 1)
    diff_otu <- append(diff_otu, diff_child)
}
diff_otu == diff_leaf

#############################################################################
