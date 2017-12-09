import matplotlib.pyplot as plt
import pickle as pkl
from PlotFunctions import *
from AAcodes import *

df = pkl.load(open('dist_hbonds.pkl','rb'))


counts = {} # dict where key = bin, value = counts
for i in [-5,-4,-3,-2,-1,1,2,3,4,5]: # fill in dict with keys
    counts[i] = 0

# store AA identity info to plot heatmap
ifgvdmdict = {} # key = bin (same as above),val = dict where (key = ifg or vdm, val = list of those)
for i in [-5,-4,-3,-2,-1,1,2,3,4,5]: # fill in dict with keys
    ifgvdmdict['ifg'] = []
    ifgvdmdict['vdm'] = []

def inc_counts_dict(counts, num):
    '''Takes a num (bin) and increases the count for that in the counts_dict
    Inputs: counts dictionary, num (bin)'''
    counts[int(num)] += 1
    return counts


def inc_heatmap_dict(ifgvdmdict, num, ifg, vdm):
    '''Takes a num (bin) and goes into that bin for the dict, and adds the 
    ifg and vdm identities (from appropriate csv file!) to the ifg and vdm lists
    Inputs: ifgvdmdict, num (bin, int), ifg (str), vdm (str)'''
    ifgvdmdict[num]['ifg'].append(ifg)
    ifgvdmdict[num]['vdm'].append(vdm)
    return ifgvdmdict

def ifg_vdm_name(row, num):
    '''Determines ifg and vdm identities from a row of a df (ex: of dist hbond vdms)
    containing columns from vdm_pdb_info csv.  
    Uses num (bin) relative resnum to get the correct amino acid from the vdm string
    Inputs: row, num (bin)'''
    ifg = row['resname_ifg']
    ifg = three_to_one[ifg]
    vdm_string = row['sequence_vdm']
    rel_resnum = num  # rel resnum for h bond, aka bin # 
    # 


def correlate_relresnum_sequence(resnums):
    '''Get index of a relresnum and then pull the AA from 
    the sequence string'''
    resnums_list = []
    if resnums.startswith('-'):
        while resnums.startswith('-'):
            resnums_list.append(int(resnums[:2]))
            resnums = resnums[2:]
    while len(resnums) > 0:
        resnums_list.append(resnums[0])
        resnums = resnums[1:]
    return resnums_list    

for ix, row in df.iterrows():
    resnums = row['rel_resnums_y']
    if len(resnums) > 1: # otherwise it's just '0'
        if resnums.startswith('-'):
            while resnums.startswith('-'):
                num = resnums[:2]
                counts = inc_counts_dict(counts, num)
                resnums = resnums[2:]
        if resnums.startswith('0'):
            while resnums.startswith('0'):
                resnums = resnums[1:] # ignore 0's
        while len(resnums) > 0:
            num = resnums[0]
            counts = inc_counts_dict(counts, num)
            resnums = resnums[1:]
            
## plot
#plt.bar(list(counts.keys()), counts.values())
#plt.title('Sequence proximity preference for H-bonds')
#plt.xlabel('relative seq #')
#plt.ylabel('counts')
#
#plt.show()
#
#


