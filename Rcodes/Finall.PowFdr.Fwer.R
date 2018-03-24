#######################################################################################
#####   This R document is to do the compare between different methods  ###############
#####   Date : 22/3/2018                                                ###############
#####   Author : ZHOU Chao                                              ###############
#####   Description :                                                   ###############
#####      Here we compare five methods on different simulation data by ###############
#####   Power and FDR                                                   ###############
#####-------------------------------------------------------------------###############
#####   Methods include :                                               ###############
#####         1. t.test on relative data ;                              ###############
#####         2. wilcox.test on count data ;                            ###############
#####         3. metagenomeSeq::ZIG on count data which will be         ###############
#####                             normalized by quanitile methods;      ###############
#####         4. ancom.R::ANCOM which use the log count ratio data      ###############
#####                               (zero will added 1);                ###############
#####         5. Tree.Ratio which based on the ratio of tree structure. ###############
#####-------------------------------------------------------------------###############
#####   Simulation:                                                     ###############
#####         1. Multinomial(Binomial) Tree distribution (BIT);         ###############
#####         2. Dirichlet Multinomial distribution (DTM)               ###############
#####         3. Logistical Normal Multinomial distribution (LNM)       ###############
#####         4. Poission distribution on each OTU introduce by ANCOM;  ###############
#####   This four distribution will fixed to one tree with OTU number   ###############
#####   p = 50, 100, 200, 500;  Each p will replicate 10 times with     ###############
#####   different tree structure, the differential proporation ranges   ###############
#####   5%-30%. Sample size case:control = 50:50 and sample depth is    ###############
#####   is 10p-1000p. Each case replicate 50 times.                     ###############
#####                  .....                                            ###############
#######################################################################################

rm(list = ls())
library(ancom.R)
library(metagenomeSeq)
source("~/Rcode/Tree.ratio.R")
source("~/Rcode/Taxa.index.R")
source("~/Rcode/Zig.pv.R")
source("~Rcode/pow.fdr.R")

sim.types <- c("BIT", "DTM", "LNM", "ANCOM")
otu.numbers <- c(50, 100, 200, 500)
tree.numbers <- 10

start <- Sys.time()
for(o in 1:length(otu.numbers)){
  p <- otu.numbers[o]
  R <- ifelse(p == 500, 30, 50)       # for big p(500) only replicate 30 times
  path.p <- paste0("~/p=", p, "/")

  #PowFdr.Sumtree <- matrix(0, 2*length(sim.types), 5)
  for(t in 1:tree.numbers){
    load(paste0(path.p, "tree_info/tree", t, ".rda"))
    dif.otu <- phy_tree$diff_otu
    taxa.index <- Taxa_index(p, phy_tree)

    PowFdr.s <- matrix(, 2*length(sim.types), 5)
    for(s in 1:length(sim.types)){
      type <- sim.types[s]

      Power.r <- Fdr.r <- matrix(, R, 5)
      for(r in 1:R){

        #------------data prepare
        path.data <- paste0(path.p, type, "/tree", t, "/otu_table.", r, ".txt")
        otu.tab <- read.table(path.data)
        colnames(otu.tab) <- as.character(1:p)
        taxa.tab <- as.matrix(otu.tab) %*% taxa.index
        all.tab <- cbind(as.matrix(otu.tab), taxa.tab)
        group <- c(rep(1, nrow(otu.tab)/2), rep(2, nrow(otu.tab)/2))

        #------------t.test
        t.rel.pv <- apply(log((otu.tab + 1e-20)/rowSums(otu.tab + 1e-20)), 2,
                          function(x) return(ifelse(length(unique(x)) == 1, 1,
                          t.test(x[group == 1], x[group == 2])$p.value)))
        t.rel.pv.adj <- p.adjust(t.rel.pv, method = "fdr", n = p)
        t.detected <- which(t.rel.pv.adj <= 0.05)

        Power.r[r, 1] <- length(intersect(dif.otu, t.detected))/length(dif.otu)
        Fdr.r[r, 1] <- ifelse(length(t.detected) == 0, 0,
                              length(setdiff(t.detected, dif.otu))/length(t.detected))

        #-------------wilcox.test
        w.count.pv <- apply(otu.tab, 2, function(x) return(wilcox.test(x[group == 1],
                                                                      x[group == 2])$p.value))
        w.count.pv.adj <- p.adjust(w.count.pv, method = "fdr", n = p)
        w.detected <- which(w.count.pv.adj <= 0.05)
        Power.r[r, 2] <- length(intersect(dif.otu, w.detected))/length(dif.otu)
        Fdr.r[r, 2] <- ifelse(length(w.detected) == 0, 0,
                              length(setdiff(w.detected, dif.otu))/length(w.detected))

        #-----------------zig
        zig.pv.adj <- Zig_pv(otu.tab, group = group)
        zig.detected <- which(zig.pv.adj <= 0.05)
        Power.r[r, 3] <- length(intersect(dif.otu, zig.detected))/length(dif.otu)
        Fdr.r[r, 3] <- ifelse(length(zig.detected) == 0, 0,
                              length(setdiff(zig.detected, dif.otu))/length(zig.detected))

        #-----------------ancom
        #ancom.res <- ANCOM(cbind(otu.tab, group), sig = 0.05, multcorr = 2)
        #ancom.deteced <- ancom.res$detected
        #dif.otux = paste0("X", dif.otu)
        #Power.r[r, 4] <- length(intersect(dif.otux, ancom.deteced))/length(dif.otux)
        #Fdr.r[r, 4] <- ifelse(length(ancom.deteced) == 0, 0,
        #                    length(setdiff(ancom.deteced, dif.otux))/length(ancom.deteced))

        #----------------tree.ratio
        tree.ratio <- Tree_ratio(p = p, tree = phy_tree, all_table = all.tab,
                                       group = group, taxa_index = taxa.index)
        tree.ratio.detected <- Tree_ratio_back(p = p, tree_ratio = tree.ratio,
                               otutab = otu.tab, taxa_index = taxa.index, group = group)
        Power.r[r, 5] <- length(intersect(dif.otu, tree.ratio.detected))/length(dif.otu)
        Fdr.r[r, 5] <- ifelse(length(tree.ratio.detected) == 0, 0,
                            length(setdiff(tree.ratio.detected, dif.otu))/length(tree.ratio.detected))

      }

      ancom <- read.table(paste0(path.p, type, "/tree", t, "/ancom.csv"), sep = ",")
      ancom.tf = ifelse(ancom[seq(2, 2*R, 2), -1] == "True", 1, 0)
      ancom.pwfdr = apply(ancom.tf, 1, function(x) return(pow.fdr(1-x, dif.otu, cutf = 0.01)))
      Power.r[, 4] = ancom.pwfdr[1, ]
      Fdr.r[, 4] = ancom.pwfdr[2, ]

      write.table(Power.r, paste0(path.p, "tree", t, "_", type, "_power.txt"))
      write.table(Fdr.r, paste0(path.p, "tree", t, "_", type, "_fdr.txt"))
      PowFdr.s[2*s-1, ] <- colMeans(Power.r, na.rm = T)
      PowFdr.s[2*s,   ] <- colMeans(Fdr.r, na.rm = T)
    }
    colnames(PowFdr.s) <- colnames(PowFdr.Sumtree) <- c("t.test", "Wilcox.test",
                                                       "Zig", "ANCOM", "Tree.Ratio")
    rownames(PowFdr.s) <- rownames(PowFdr.Sumtree) <- c("BIT.Pwr", "BIT.Fdr", "DTM.Pwr",
                         "DTM.Fdr", "LNM.Pwr", "LNM.Fdr", "ANCOM.Pwr", "ANCOM.Fdr")
    write.table(PowFdr.s, paste0(path.p, "MeanPowerFdr_", "tree", t, ".txt"))
    PowFdr.Sumtree <- PowFdr.Sumtree + PowFdr.s  #
  }
  PowFdr.MeanTree <- PowFdr.Sumtree/tree.numbers
  write.table(PowFdr.MeanTree, paste0(path.p, "POW.FDR.MeanTree.txt"))
}

Sys.time() - start

#--------------------------------------------------------------------------------------------#

  ###########################################################################################
 ##    Next, we compare this methods by FWER(famliy wise error rate) which calculate by     ##
##  permutate the group label randomly L = 10 times for each data, this will be match to the ##
##  null hypothesis, so once reject one show the method happen to type I error, take average ##
 ##  of this value by each permutation and data get the FWER.                               ##
 ##                                                                                         ##
 #############################################################################################

L <- 100
R <- 10   # choose 10 data to do this

start <- Sys.time()
for(o in 1:length(otu.numbers)){
  p <- otu.numbers[o]
  path.p <- paste0("~/p=", p, "/")

  FWER.Sumtree <- matrix(0, length(sim.types), 5)
  for(t in 1:tree.numbers){
    load(paste0(path.p, "tree_info/tree", t, ".rda"))
    taxa.index <- Taxa_index(p, phy_tree)

    FWER.type <- matrix(, length(sim.types), 5)
    for(s in 1:length(sim.types)){
      type <- sim.types[s]
      FWER.r <- matrix(NA, nrow = L*R, ncol = 5)
      colnames(FWER.r) <- c("t.test", "wilcox.test", "Zig", "ANCOM", "Tree.Ratio")
      group.permt <- matrix(1, nrow = L*R, ncol = ifelse(p == 500, 60, 100))

      for(r in 1:R){
        #------------data prepare
        path.data <- paste0(path.p, type, "/tree", t, "/otu_table.", r, ".txt")
        otu.tab <- read.table(path.data)
        colnames(otu.tab) <- as.character(1:p)
        taxa.tab <- as.matrix(otu.tab) %*% taxa.index
        all.tab <- cbind(as.matrix(otu.tab), taxa.tab)

        for(l in 1:L){
          group = rep(1, nrow(otu.tab))
          set.seed(r*s*t*o*l)
          group[sample(1:nrow(otu.tab), nrow(otu.tab)/2)] <- 2
          group.permt[(r-1)*L + l, ] <- group

          #------------t.test
          t.rel.pv <- apply(log((otu.tab + 1e-20)/rowSums(otu.tab + 1e-20)), 2,
                          function(x) return(ifelse(length(unique(x)) == 1, 1,
                                                    t.test(x[group == 1], x[group == 2])$p.value)))
          t.rel.pv.adj <- p.adjust(t.rel.pv, method = "fdr", n = p)
          FWER.r[(r-1)*L + l , 1] <- ifelse(sum(t.rel.pv.adj <= 0.05) != 0, 1, 0)

          #-------------wilcox.test
          w.count.pv <- apply(otu.tab, 2, function(x) return(wilcox.test(x[group == 1],
                                                                       x[group == 2])$p.value))
          w.count.pv.adj <- p.adjust(w.count.pv, method = "fdr", n = p)
          FWER.r[(r-1)*L + l, 2] <- ifelse(sum(w.count.pv.adj <= 0.05) != 0, 1, 0)

          #-----------------zig
          zig.pv.adj <- Zig_pv(otu.tab, group = group)
          zig.detected <- which(zig.pv.adj <= 0.05)
          FWER.r[(r-1)*L + l, 3] <- ifelse(length(zig.detected) != 0, 1, 0)

          #-----------------ancom
          ancom.res <- ANCOM(cbind(otu.tab, group), sig = 0.05, multcorr = 2)
          FWER.r[(r-1)*L + l, 4] <- ifelse(sum(ancom.res$detected != "No significant OTUs detected") == 0
                                           , 1, 0)

          #----------------tree.ratio
          tree.ratio <- Tree_ratio(p = p, tree = phy_tree, all_table = all.tab,
                                 group = group, taxa_index = taxa.index)
          tree.ratio.detected <- Tree_ratio_back(p = p, tree_ratio = tree.ratio,
                                               otutab = otu.tab, taxa_index = taxa.index, group = group)
          FWER.r[(r-1)*L + l, 5] <- ifelse(length(tree.ratio.detected) != 0, 1, 0)
        }

      }

      # # use ancomP in python which is more faster than in R
      # ancom.permt.result <- read.table(paste0(path.p, type, "/tree", t, "/ancom.permt.csv", sep = ","))
      # FWER.r[, 4] <- apply(ancom.permt.result[seq(3, (3*L), 3), -1], 1, function(x)
      #                                       return(ifelse(sum(x == "True") != 0, 1, 0)))
      #
      write.table(FWER.r, paste0(path.p, "tree", t, "_", type, "_fwer.txt"))
      write.table(group.permt, paste0(path.p, type, "/tree", t, "/group.permt.txt"),
                  col.names = F, row.names = F)

      colMeans(FWER.r, na.rm = T) -> FWER.type[s, ]
    }
    colnames(FWER.type) <- colnames(FWER.Sumtree) <- c("t.test", "Wilcox.test",
                                                        "Zig", "ANCOM", "Tree.Ratio")
    rownames(FWER.type) <- rownames(FWER.Sumtree) <- c("BIT", "DTM", "LNM", "ANCOM")

    write.table(FWER.type, paste0(path.p, "MeanFWER_", "tree", t, ".txt"))
    FWER.Sumtree <- FWER.Sumtree + FWER.type
  }

  FWER.MeanTree <- FWER.Sumtree/tree.numbers
  write.table(FWER.MeanTree, paste0(path.p, "FWER.MeanTree.txt"))
}

Sys.time() - start

# save.image(file = "throat.RData")


# permutate the group label of each data
start <- Sys.time()
for(o in 1:length(otu.numbers)){
  p <- otu.numbers[o]
  path.p <- paste0("~/p=", p, "/")
  for(t in 1:tree.numbers){
    for(s in 1:length(sim.types)){
      type <- sim.types[s]
      N = ifelse(p == 500, 60, 100)
      group.permt <- matrix(1, nrow = L*R, ncol = N)

      for(r in 1:R){
        for(l in 1:L){
          group = rep(1, N)
          set.seed(r*s*t*o*l)
          group[sample(1:N, N/2)] <- 2
          group.permt[(r-1)*L + l, ] <- group
        }
      }
      write.table(group.permt, paste0(path.p, type, "/tree", t, "/group.permt.txt"),
                  col.names = F, row.names = F)

    }
  }
}
Sys.time() - start
