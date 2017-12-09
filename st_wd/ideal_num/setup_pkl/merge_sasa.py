import sys                              
import pandas as pd
import pickle as pkl
from sys import argv

script, ifg, csv_dir = argv

df = pkl.load(open('noskip_%s_sasa.pkl'%ifg,'rb'))
df = df[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'vdM_sasa_dssp_full_residue', 'vdM_sasa_CB_3A_probe']]

# add atom names and hbond info
atoms = pd.read_csv(csv_dir+'%s_vdm_pdb_info.csv'%ifg)
atoms = atoms[['iFG_count', 'vdM_count', 'atom_names']]
merged = pd.merge(df,atoms, on=['iFG_count', 'vdM_count'])

# add vdm hbond info
def add_hbond(merged, hb):
    ''' hb can be 'hbond' or 'ca_hbond' '''
    vdmhbond = pd.read_csv(csv_dir+'%s_ifg_%s_vdm.csv'%(ifg, hb))[['iFG_count', 'vdM_count', 'rel_resnums', 'number_hbonds']] # number_hbonds is only hbonds with vdm
    if hb == 'ca_hbond':
        vdmhbond = vdmhbond.rename(columns={'number_hbonds':'ca_hbonds'})
        vdmhbond = vdmhbond.loc[vdmhbond['rel_resnums'].isin(['0', '00', '000', '0000'])][['iFG_count', 'vdM_count', 'ca_hbonds']]
    else:
        vdmhbond = vdmhbond.loc[vdmhbond['rel_resnums'].isin(['0', '00', '000', '0000'])][['iFG_count', 'vdM_count', 'number_hbonds']]
    merged = pd.merge(merged, vdmhbond, on=['iFG_count','vdM_count'], how='left')
    return merged

merged = add_hbond(merged, 'hbond')
merged = add_hbond(merged, 'ca_hbond')

pkl.dump(merged, open('noskip_%s_contacts.pkl' %ifg, 'wb'))
