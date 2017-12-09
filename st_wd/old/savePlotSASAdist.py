# commandline argument is fs, dssp, 3A, 4A, or 5A
import numpy as np
import pickle as pkl
from pprint import pprint
from sys import argv
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

script,m = argv

nrpdb_dict=pkl.load(open('nrPDB_sasa.pkl','rb'))

MaxASA = {'ALA': 129, 'ARG': 274, 'ASN': 195, 'ASP':193, 'CYS':167, 'GLU':223, 'GLN':225, 'GLY':104,
          'HIS':224, 'ILE':197, 'LEU':201, 'LYS':236, 'MET':224, 'PHE':240, 'PRO':159, 
          'SER':155, 'THR':172, 'TRP':285, 'TYR':263, 'VAL':174} 

f, axarr = plt.subplots(4,5,figsize=(12,6.5), sharey=True, sharex=True)
fig_row = 0
fig_col = 0

for aa, resdict in nrpdb_dict.items():
    x = resdict[m]
    x = [float(z) for z in x ]
    x = np.hstack(x)
    #x = x/MaxASA[aa]
    x = sorted(x)


    axarr[fig_row, fig_col].hist(x, bins=10)
    axarr[fig_row, fig_col].set_title(aa)

    if fig_col < 4:
        fig_col += 1
    else:
        fig_col = 0
        fig_row += 1

plt.xlim(0,250)
f.text(0.5, 0.95, 'nrPDB sasa cb probe size %s' %m, ha='center')
f.text(0.5, 0.04, 'sasa sq. angstroms', ha='center')
f.text(0.04, 0.5, 'counts', va='center', rotation='vertical')
#plt.tight_layout()
plt.show()
###
