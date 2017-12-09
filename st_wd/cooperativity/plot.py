import pickle as pkl, pandas as pd, seaborn as sns, matplotlib.pyplot as plt
import numpy as np
from AAcodes import *
from sys import argv
import scipy

script, correlation_pkl = argv
corr = pkl.load(open(correlation_pkl,'rb'))
scores,obs_exp = corr
ifg = correlation_pkl.split('_')[0]
# make heatmap
aas = sorted(list(one_to_three.keys()))
heatmap = pd.DataFrame(index=aas, columns=aas)
heatmap = heatmap.replace(np.NaN, 0)
heatmap = heatmap.astype(float)

# generate mask for lower triangle
mask = np.zeros_like(heatmap, dtype=np.bool)
mask[np.triu_indices_from(mask,k=1)] = True

for ix1, AA1 in enumerate(aas): # aa1 = aai
    for ix2, AA2 in enumerate(aas):
        aa1, aa2 = one_to_three[AA1], one_to_three[AA2]
        try:
            heatmap[AA1][AA2] = np.round(scores[(aa1,aa2)],2)*10
        except:
            heatmap[AA1][AA2] = np.round(scores[(aa2,aa1)],2)*10
            
        
        mx = np.ma.masked_array(heatmap, mask=mask)
vmax,vmin,cen = mx.max(), mx.min(), mx.mean()
#print(vmax,vmin,cen)

sns.heatmap(heatmap, annot=True,mask=mask)
#sns.heatmap(heatmap, annot=True,mask=mask, cmap="YlGnBu", center=0,vmin=vmin,
#        square=True, linewidths=.5, cbar_kws={"shrink": .5})
#sns.heatmap(heatmap, annot=True,mask=mask, cmap="YlGnBu", vmax=vmax,center=cen,
#        vmin=vmin,    square=True, linewidths=.5, cbar_kws={"shrink": .5})
ticks = np.arange(20) + 0.5
plt.yticks(ticks,aas[::],rotation=0) 
plt.xticks(ticks,aas,rotation=0) 
plt.title('%s cooperativity score * 10'%ifg)
plt.show()
