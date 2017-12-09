import seaborn as sns
import matplotlib.pyplot as plt
from sys import argv
import pandas as pd
from PhiPsiFunctions import *

script, carboxamide_df_path, db_df_path = argv

carboxamide = log_norm_df(carboxamide_df_path)
db = log_norm_df(db_df_path)

df = carboxamide/db

# x axis = col (phi), y axis = index (psi)
sns.heatmap(df, yticklabels = False, xticklabels=False)
plt.title('difference between normalized log(counts) of carboxamide vdm and database AA phi/psi')
plt.show()
