from pprint import pprint
import matplotlib.pyplot as plt
import pickle as pkl
from PlotFunctions import *
from AAcodes import *
import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.analysis import Analysis


df = pkl.load(open('dist_hbonds.pkl','rb'))
ind = Analysis.repeat_indices()
df = pd.DataFrame(df, index=ind)

#df = df[['pdb','resname_ifg', 'resname_vdm', 'sequence_vdm', 'sec_struct_dssp_vdm']]
x = Analysis.remove_repeat_proteins(df)

#
#
#counts = {} # dict where key = bin, value = counts
#for i in [-5,-4,-3,-2,-1,1,2,3,4,5]: # fill in dict with keys
#    counts[i] = 0
#
## store AA identity info to plot heatmap
#heatmapdict = {} # key = bin (same as above),val = dict where (key = vdm0 or vdmi, val = list of those)
#for i in [-5,-4,-3,-2,-1,1,2,3,4,5]: # fill in dict with keys
#    heatmapdict[i] = {}
#    heatmapdict[i]['vdm0'] = []
#    heatmapdict[i]['vdmi'] = []
#    heatmapdict[i]['index'] = []
#
#for ix, row in df.iterrows():
#    hb_resnums = parse_resnums(row['rel_resnums_y']) # for resnums in vdmfrag ifg interacts with 
#    for num in hb_resnums:
#        if num != 0:
#            try:
#                # add to counts dict for histogram plot
#                counts = inc_counts_dict(counts, num)
#                # add to heatmapdict for heatmap
#                vdm0, vdmi = vdm0_vdmi_names(row, num)
#                if vdm0 == 'L' and vdmi == 'L' and num == 5:
#                    print(ix)
#                heatmapdict = inc_heatmap_dict(heatmapdict, num, vdm0, vdmi)
#            except:
#                pass
#
#x = (heatmapdict[5])
##for i in range(len(x['vdm0'])):
##
##    if x['vdm0'][i]=='L' and x['vdmi'][i]=='L':
##        print(i)
##
#
#
#
### plot
##plt.bar(list(counts.keys()), counts.values())
##plt.title('Sequence proximity preference for H-bonds')
##plt.xlabel('relative seq #')
##plt.ylabel('counts')
##
##plt.show()
##
##
##
##
