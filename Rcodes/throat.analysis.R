###-----------------Ratio test incroporate the phylogeney information apply to throat data----------###
##
##  This code main to compare different methods including (t.test(relative), zig(metagenomeSeq),
## ANCOM(ancom.R), tree.ratio(Tree.Ratio) and NonTree.Ratio (Same methods just like ANCOM) on the throat
## data from Chong WU (https://github.com/ChongWu-Biostat/MiSPU/tree/master/data)

# library nacessary packages and load the data
rm(list = ls())
library(MiSPU)
library(ancom.R)
library(metagenomeSeq)
source("~/Rcode/Tree.ratio.R")
source("~/Rcode/Taxa.index.R")
source("~/Rcode/NonTree.test.R")
source("~/Rcode/Zig.pv.adj.R")
data("throat.otu.tab")
data("throat.meta")
data("throat.tree")

###------------------------------------------- pre-process ----------------------------------------###
## this is a binary tree which is suitable for our methods
## Because the OTU numbers is a little big, so here take 0.1 as the significant level
threshold <- 0.1
otu.number <- ncol(throat.otu.tab)
is.binary.tree(throat.tree)
colnames(throat.otu.tab) <- as.character(1:otu.number)
throat.group <- which(throat.meta[,3] == "Smoker")


###-------------------- t.test of the relative otu.table data for differential depth---------------###
throat.log.rel.data = log((throat.otu.tab + 1e-20)/rowSums(throat.otu.tab + 1e-20))
throat.t.rel.pvalue <- apply(throat.log.rel.data, 2, function(x)
  return(ifelse(length(unique(x)) == 1, 1, t.test(x[throat.group], x[-throat.group])$p.value)))
throat.t.rel.pvalue.adj <- p.adjust(throat.t.rel.pvalue, method = "fdr", n = otu.number)
throat.t.rel.detected <-  names(which(throat.t.rel.pvalue.adj <= threshold))


###-------------------------- ZIG which is use the quantile normalization -------------------------###
rownames(throat.otu.tab) <- as.character(1:nrow(throat.otu.tab))
throat.zig.group <- AnnotatedDataFrame(data.frame(throat.meta$SmokingStatus))
throat.zig.otu.tab <- newMRexperiment(t(throat.otu.tab), phenoData = throat.zig.group)
throat.zig.otu.norm <- cumNorm(throat.zig.otu.tab, p = cumNormStatFast(throat.zig.otu.tab))
throat.zig.model <- model.matrix(~ 1 + throat.meta$SmokingStatus,
                                 data = pData(throat.zig.otu.norm))
throat.zig.fit <- fitFeatureModel(throat.zig.otu.norm, throat.zig.model)
throat.zig.result <- MRcoefs(throat.zig.fit, number = otu.number)
throat.zig.detected <- rownames(throat.zig.result[which(throat.zig.result[, 4] <= threshold ), ])


###----------------- ANCOM which take logarithm ratio transform for each OTU ---------------------###
## Have to infrom that use ancom.R::ANCOM which in this data will consume several minutes, just be
## PATIENT!!! (Good Luck to you for more OTU numbers data use this function), You can also try to
## use ancomP in python which is much more faster.
## Reference :1) https://www.niehs.nih.gov/research/resources/software/biostatistics/ancom/index.cfm
##            2) https://github.com/alk224/akutils-v1.2/blob/master/akutils_resources/R-instructions_ancom.r
start = Sys.time()
throat.ancom.otu.tab <- cbind(throat.otu.tab, throat.meta$SmokingStatus)
throat.ancom.result <- ANCOM(throat.ancom.otu.tab, sig = 0.05, multcorr = 2)
Sys.time() - start
throat.ancom.detected <- throat.ancom.result$detected

# here you can try to use the python code to do this
# Reference: https://github.com/mortonjt/ancomP
write.table(throat.otu.tab, "/Users/zhouchao/Documents/Tree/OTU/throat.otu.tab.txt",
            row.names = F, col.names = F)
write.table(t(throat.meta$SmokingStatus), "/Users/zhouchao/Documents/Tree/OTU/throat.group.txt",
            row.names = F, col.names = F, sep = ",")
throat.ancom.py.result <- read.table("/Users/zhouchao/Documents/Tree/throat/ancom.result.csv",
                                     sep = ",", head = T)
throat.ancom.py.delected <- which(throat.ancom.py.result[, 3] == "True")


###--------------- NonTree take the reject proportion to do decesion of p-1 otu ------------------###
## Although ancom.R::ANCOM use the parallel technolongy in R but it is still very slow, nontree ratio
## test is just the same with it. (this spend about 2.4 min)
throat.nontree.result <- NonTree_test(throat.otu.tab, p = otu.number, group = throat.group, log = T)
throat.nontree.detected <- which(throat.nontree.result[, 3] > 0.5)

###---------------------------- Tree Ratio test with correction ----------------------------------###
throat.taxa.index <- Taxa_index(p = otu.number, phy_tree = throat.tree)
throat.taxa.table <- as.matrix(throat.otu.tab)%*%throat.taxa.index
throat.all.table <- cbind(as.matrix(throat.otu.tab), throat.taxa.table)
throat.tree.ratio <- Tree_ratio(p = otu.number, tree = throat.tree, all_table = throat.all.table,
                               group = throat.group, taxa_index = throat.taxa.index)
throat.taxa.diff <- throat.tree.ratio$taxa.dif
throat.tree.back <- Tree_ratio_back(p = otu.number, tree_ratio = throat.tree.ratio,
                                   otutab = throat.otu.tab, group = throat.group)
throat.tree.back.adj <- p.adjust(throat.tree.back, method = "fdr", n = length(throat.tree.back))
throat.tree.ratio.detected <- names(which(throat.tree.back.adj < 0.1))


###------------finally we can get the differential OTU detected by this methods--------------------###
## we can plot venn plot which show that the tree.ratio test can find the most OTU which content the
## OTU detected by ANCOM while t.test and zig have no detected or just 1. However, we can use the FWER
## (Family-Wise Error Rate) to compare the type I error of this methods. So, we permutate the label of
## the Smoking Status 1000 times to calculate it as follow.

L <- 100
throat.FWER <- matrix(NA, nrow = L, ncol = 5)
colnames(throat.FWER) <- c("t.test", "Zig", "NonTree", "ANCOM", "Tree.Ratio")
throat.grouping.permt.label <- matrix(1, nrow = L, ncol = nrow(throat.meta))

start.permt = Sys.time()
for(i in 1:L){
  set.seed(i)
  throat.grouping.permt <- sample(1:nrow(throat.meta), 28)
  throat.grouping.permt.label[i, throat.grouping.permt] <- 2

  # t.test
  throat.t.rel.pv.permt <- apply(throat.log.rel.data, 2, function(x)
    return(ifelse(length(unique(x)) == 1, 1, t.test(x[throat.grouping.permt], x[-throat.grouping.permt])$p.value)))
  throat.t.rel.p.adj.permt <- p.adjust(throat.t.rel.pv.permt, method = "fdr", n = otu.number)
  throat.FWER[i, 1] <- sum(throat.t.rel.p.adj.permt <= 0.05)

  # zig, here functional the zig process
  throat.zig.p.adj.permt <- Zig.pv.adj(throat.otu.tab, group = throat.grouping.permt.label[i, ])
  throat.FWER[i, 2] <- sum(throat.zig.p.adj.permt[-which(is.na(throat.zig.p.adj.permt))] <= 0.05)

  # nontree
  throat.nontree.rp.permt <- NonTree.test(throat.otu.tab, p = otu.number,
                                              group = throat.grouping.permt.label[i, ], log = T)[, 3]
  throat.FWER[i, 3] <- sum(throat.nontree.rp.permt >= 0.5)

  # tree.ratio
  throat.tree.ratio.permt <- Tree.ratio(p = otu.number, tree = throat.tree, all_table = throat.all.table,
                                  group = throat.grouping.permt, taxa_index = throat.taxa.index)
  throat.tree.ratio.pv.permt <- Tree.ratio.back(p = otu.number, tree_ratio = throat.tree.ratio.permt,
                                      otutab = throat.otu.tab, taxa_index = throat.taxa.index, group = throat.grouping.permt)
  throat.FWER[i, 5] <- sum(throat.tree.ratio.pv.permt <= 0.05)

}
Sys.time() - start

# use ancomP in python which is more faster than in R
write.table(throat.grouping.permt.label, "throat.grouping.permt.label.txt", row.names = F, col.names = F)
ancom.permt.result <- read.table("ancom.permt.result.csv", sep = ",")
throat.FWER[, 4] <- apply(ancom.permt.result[seq(3, (3*L), 3), -1], 1, function(x)
                                                      return(ifelse(sum(x == "True") != 0, 1, 0)))

write.table(throat.FWER, "throat.FWER.txt", col.names = T, row.names = F)

throat.FWER[throat.FWER != 0] <- 1
print(colMeans(throat.FWER))

setwd("/Users/zhouchao/Documents/Tree/throat"); getwd()
save.image(file = "throat.RData")







