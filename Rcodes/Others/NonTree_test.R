#Just like ANCOM, but take the secondary minimum and maximum pvalue or reject proportion to do decision.
#Here call this function NonTree test
##ancom + sed-min -max

NonTree.test = function(otu_table, p, group = c(1:2*N), log = F){
    
  ancom_table <- matrix(NA, p, p-1)
  otu_table <- as.matrix(otu_table)
  group1 <- which(group == unique(group)[1])
  
  for(i in 1:p){
    if(log == F){
      temp_table <- (otu_table[, i]+1e-30)/(otu_table[, -i]+1e-30)
    }
    
    else{
      temp_table <- log(otu_table[, i] + 1e-30) - log(otu_table[, -i] + 1e-30)
    }
    ancom_table[i, ] <- apply(temp_table, 2, function(x) return(ifelse(length(unique(x)) == 1, 1,
                                                  t.test(x[group1], x[-group1])$p.value)))
  }
  ancom_m2min <- apply(ancom_table, 1, function(x) return(sort(na.omit(x))[2]))
  ancom_m2min.adj <- p.adjust(ancom_m2min, method = "fdr", n = p)
  
  ancom_m2max <- apply(ancom_table, 1, function(x) return(sort(na.omit(x))[length(na.omit(x))-2]))
  ancom_m2max.adj <- p.adjust(ancom_m2max, method = "fdr", n = p)
  
  ancom_rejectratio <- apply(ancom_table, 1, function(x) return(sum(na.omit(x)<=0.05)/(p-1)))
  result 《- data.frame( NonTree_SecMax = ancom_m2max.adj, NonTree_SecMin = ancom_m2min.adj,
                       NonTree_Reject = ancom_rejectratio)
  return(result)
}
# = NonTree_test(otu_table, p , index = c(1:N))
