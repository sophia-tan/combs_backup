# generate lookup table for counts of AAi in combs output
import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis
from sys import argv
import pickle as pkl, pandas as pd

script, csv_dir, ifg = argv

an = analysis.Analyze(csv_dir)
df = an.get_distant_vdms(seq_distance=10)
print(len(df), 'number of vdMs w/ skipping #10')
df = analysis.Analysis.remove_repeat_proteins(df)

# add info about ifg-vdm dist from contacts csv
contactsdf = an.ifg_contact_vdm
merged = pd.merge(df,contactsdf, on=['iFG_count', 'vdM_count'])
merged = merged[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'atom_names_vdm', 'sequence_vdm', 'sec_struct_dssp_vdm', 'dist_info']]

print(len(merged), 'number of vdMs w/ repeats removed')
aai_df = analysis.EnergyTerms.make_freqaai_df(merged)
pkl.dump(aai_df, open('AAi_freq_combed_%s.pkl' % ifg, 'wb'))
