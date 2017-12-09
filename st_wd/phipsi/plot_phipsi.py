import pickle as pkl
import seaborn as sns
import matplotlib.pyplot as plt
from sys import argv
import numpy as np
import pandas as pd
from PhiPsiFunctions import *

script, df_path = argv

df = log_norm_df(df_path)

# x axis = col (phi), y axis = index (psi)
sns.heatmap(df, yticklabels = False, xticklabels=False)
plt.title('normalized log(counts) of database AA phi/psi')
#plt.title('normalized log(counts) of carboxamide vdm phi/psi')
plt.show()
