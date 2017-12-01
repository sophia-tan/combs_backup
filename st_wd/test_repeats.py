import pickle as pkl
import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis

df = pkl.load(open('armadillo.pkl','rb'))
removed = pkl.load(open('removed.pkl','rb'))

df = analysis.Analysis.remove_repeat_proteins(df)

removed = removed.index
df = df.index

