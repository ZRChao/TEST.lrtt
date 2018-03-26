## Tree.Ratio
Differential Abundance Analysis for Microbiome data Incorporating Phylogeny

## Installation

```R
# Install the released version from CRAN
install.packages("LTRR")

# Install the development version from GitHub
devtools::install_github("ZRChao/LTRR")
```
## Contents

### Simulation
* BIT.Sim.R for Multinomial(Binomial) Tree distribution
* DTM.Sim.R for Dirichlet Multinomial Tree distribution
* LNM.Sim.R for Logistical Normal Multinomial distribution 
* ANCOM.Sim.R for Poission distribution (parameters set follow as ANCOM paper)

### Tree relate function
* Tree.Sim.R : Tree simulation
* Taxa.index.R : relationship between leafs and internal nodes
* Prob.mult.R :calculate the probability of each leafs by multiple each probability along the branch
* Diff.otu.R : judge the probability of each node in two group equally or not
* Tree.ratio.R : Tree ratio test based on the tree structure
* 

### Other function
* pow.fdr.R : calculate power and fdr with p.value 
* Zig.pv.adj.R : metagenomeSeq::fitFeature
* ancom.py : which is faster in python

### Require 
* ancom.R::ANCOM, metagenomeSeq::fitFeature, mvtnorm, gtools, rdirimult 

### Example
 
 You can see this in real data application of throat.R or gut.R

