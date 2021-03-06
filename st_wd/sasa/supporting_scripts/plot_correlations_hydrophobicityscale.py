import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from scipy import stats

carbox_df = pkl.load(open('scores_largeprobe_sasa_carboxylate_lookup.pkl','rb'))
db_df = pkl.load(open('scores_largeprobe_sasa_db_lookup.pkl','rb'))

Eisen = [0.620  ,-2.530  ,-0.780  ,-0.900  , 0.290  ,-0.850  ,-0.740  , 0.480  ,-0.400  , 1.380  , 1.060  ,-1.500  , 0.640  , 1.190  , 0.120  ,-0.180  ,-0.050  , 0.810  , 0.260  , 1.080 ]
AAs = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
df = pd.DataFrame(index=AAs, columns=['carboxylate', 'db', 'Eisenberg'])

for index, aa in enumerate(AAs):
    df.ix[aa,'Eisenberg'] = Eisen[index]

x = []
y = []
list_df = [carbox_df, db_df]
for aa in AAs:
    if aa != 'ARG':
        for index, each in enumerate(['carboxylate', 'db']):
            first_bin = list_df[index].ix[10, aa]
            #second_bin = list_df[index].ix[20, aa]
            #avg = np.mean([first_bin, second_bin])
            #df.ix[aa, each] = avg
            df.ix[aa,each] = first_bin

df = df.dropna(axis=0, how='any')


corr_carbox = stats.pearsonr(df['carboxylate'], df['Eisenberg']) 
corr_db = stats.pearsonr(df['db'], df['Eisenberg']) 

print('pearson coefficients for carboxylate vdms and database AAs: ')
print(corr_carbox[0], corr_db[0])

for label, x, y in zip(AAs, df['carboxylate'], df['Eisenberg']):
    plt.annotate(
        label,
        xy=(x, y), xytext=(10, 10),
        textcoords='offset points', ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.2),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
for label, x, y in zip(AAs, df['db'], df['Eisenberg']):
    plt.annotate(
        label,
        xy=(x, y), xytext=(10, -10),
        textcoords='offset points', ha='center', va='top',
        bbox=dict(boxstyle='round,pad=0.1', fc='yellow', alpha=0.2),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.scatter(df['carboxylate'],df['Eisenberg'])
plt.scatter(df['db'],df['Eisenberg'])
plt.title('cb large probe sasa vs. hydrophobicity scale; b=carboxylate vdm, o=db AA')

plt.show()
