#############################################################################
###            Zig differential analysis                                #####
### Here to more convenient to use the package metagenomeSeq we package #####
### it as one function and return the adjusted pvalue of each OTU       #####
###                                                                     #####
#############################################################################

#-----------------------------------------------------------------------#####

# need to library(metagenomeSeq)
Zig.pv.adj = function(data, group){
    
  temp_data = t(as.matrix(data))
  indz = AnnotatedDataFrame(data.frame(X = group))
  datz = newMRexperiment(temp_data, phenoData = indz) \
  qua = cumNormStatFast(datz)
  new.data = cumNorm(datz, p = qua)#

  mod = model.matrix(~ 1 + X, data = pData(new.data))
  res = fitFeatureModel(new.data, mod)
  res.pv = MRcoefs(res, number = ncol(data))
  
  zig_pv = NULL
  for(z in 1:ncol(data)){
    zig_pv[z] = res.pv[which(rownames(res.pv) == as.character(z)), 4]
  }
  return(zig_pv)
}
##zig function to two group and get the p.value.adj
# data with each row is one sample with different OTU
# the otu names is as.character 1:otu.number and finally
# can get the pvalue sorted by the name
#############################################################################
