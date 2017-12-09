from pprint import pprint
import matplotlib.pyplot as plt
import pickle as pkl
from AAcodes import *
import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs.analysis import Analysis
from combs.analysis.Analysis import *
from combs.analysis import PlotFunctions
import pandas as pd
import traceback
import matplotlib
import numpy as np
#matplotlib.style.use('ggplot')

df = pkl.load(open('dist_hbonds_no_repeats.pkl','rb'))

counts_dict = {} # key = ifgcount, v = # hbonds
for ix, row in df.iterrows():
    ifg = row['iFG_count']
    hb_resnums = PlotFunctions.parse_resnums(row['rel_resnums_y'])
    hb_resnums = [x for x in hb_resnums if x==0]
    try:
        counts_dict[ifg] += len(hb_resnums)
    except:
        counts_dict[ifg] = len(hb_resnums)
values = counts_dict.values()
values = [int(x) for x in values]
print(np.mean(values))
