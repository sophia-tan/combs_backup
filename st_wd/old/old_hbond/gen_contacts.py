## need to incorporate into combs at some point!

import sys
sys.path.append('/home/gpu/Sophia/combs/src/')
import combs
import pickle as pkl
from combs import analysis
import numpy as np
import pandas as pd

an = analysis.Analyze('/home/gpu/Sophia/STcombs/20170725/csv')
# 1) get ifg_count and vdm_count of distant vdms
dist_vdms = an.get_distant_vdms(seq_distance=10)
# select dist vdms in other csv files by getting ifg and vdm counts
ifgvdmcount = zip(dist_vdms['iFG_count'], dist_vdms['vdM_count'])


## 2) use ifgvdmcount to grab dist vdms from hbdond pd
hbond_pd = an.ifg_hbond_vdm
dist_hbond_relresnums = []

for ifg, vdm in ifgvdmcount:
    dist_hbond = hbond_pd.loc[hbond_pd['iFG_count']==ifg] # check to see if same ifgcount
    dist_hbond = dist_hbond.loc[dist_hbond['vdM_count'] == vdm]
    # may or may not be in hbond_csv file! not all distant vdms make h bonds 
    if len(dist_hbond) == 1: # if that dist vdm is in hbond csv file
        for x in dist_hbond['rel_resnums']: # have to unpack it like this bc it's a series
            dist_hbond_relresnums.append(x)


with open('relresnums.pkl', 'wb') as f:
    pkl.dump(dist_hbond_relresnums, f)



#for row in dist_hbond_relresnums:
#    for rel_resnum in row:
#        if '0' in rel_resnum:
#            print(rel_resnum)



#for ifg, vdm in ifgvdmcount:
#    for row_ix, row in hbond_pd.iterrows():
#        if row['iFG_count']==ifg and row['vdM_count']==vdm:
#            print(row['rel_resnums'])


##    and grab rel_resnums column
#hbond_rows = hbond_pd.ix[indices,'rel_resnums']
#
#counts = 0
#print(len(hbond_rows))
#for i in hbond_rows:
#    if pd.isnull(i) == True:
#        counts += 1
#print(counts)
#
#
#
#
#
