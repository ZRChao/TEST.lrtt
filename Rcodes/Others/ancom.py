# This document is phthon code for ancomP.stats.ancom::ANCOM to test differential OTU
# which is much more faster than R package ancom.R::ANCOM

from ancomP.stats.ancom import ancom
import pandas as pd
import numpy as np
from skbio.stats.composition import multiplicative_replacement

p = 20

for j in range(50):
    dir1 = 'H:/Tree/tree_base/p=' str(p) + '/otu_table.' + str(j+1) +'.txt'
    data = open(dir1, 'r')
    tmp=[]
    lines = data.readlines()
    for line in lines:
        line = list(line.strip().split(' '))
        s = []
        for n in line:
            s.append(int(n))
        tmp.append(s)
    data.close()
    dat = multiplicative_replacement(tmp)
    ind = np.arange(1,p+1,1)
    sam = np.arange(1,101,1)
    table = pd.DataFrame(dat, index = sam, columns = ind)
    grouping = pd.Series(sorted([0,1]*50),index = sam)

    results = ancom(table, grouping) # default parameters
    resultsT = results.T
    resultsT.to_csv('H:/Tree/tree_base/p=' + str(p) + '/ANCOM.csv', mode = 'a',header = False)
