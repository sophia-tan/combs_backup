import pickle as pkl
from sys import argv
import numpy as np
from PhiPsiFunctions import *

script, df_path = argv
vdm_df = pkl.load(open(df_path,'rb'))

counts_df = create_phi_psi_matrix()
counts_df = inc_phi_psi_df(counts_df, vdm_df)
print(counts_df)

pkl.dump(counts_df, open('phipsi_bin_counts_df.pkl','wb'))
