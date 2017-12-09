import itertools
import numpy as np, pickle as pkl

resname_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
                   }

def correlation(df):
    ''' Returns 2 dictionaries: (1) scores, and (2) raw values for num_obs and num_exp'''
    scores = {}
    obs_exp = {}
    df = df[df.apply(get_resnames, axis=1)]
    vdm_gr = df.groupby('iFG_count')
    vdm_gr_agg_pairs = vdm_gr['resname_vdm'].agg([num_pairs, to_pairs])
    vdm_gr_agg_pairs_gte1 = vdm_gr_agg_pairs[vdm_gr_agg_pairs['num_pairs'] > 0]
    n_pairs_total = vdm_gr_agg_pairs_gte1['num_pairs'].sum()
    resns = set(df.groupby('resname_vdm').groups)
    print(n_pairs_total,'n pairs total')
    
    ### make a dictionary that gives the freq of finding a given AA in a pair ###
    ### this will be used to calculate the expectation num                    ###
    ### need to also account for if this AA is in a pair with itself!         ###
        
    def count_obs_pairs(row, a, b):
        # order doesn't matter
        if a != b:
            return row['to_pairs'].count((a, b)) + row['to_pairs'].count((b,a))
        else: 
            return row['to_pairs'].count((a, b))
            
    for res1, res2 in itertools.combinations_with_replacement(sorted(list(resns)), 2):
        n_obs_res1_res2_pairs = vdm_gr_agg_pairs_gte1.apply(count_obs_pairs,axis=1,args=(res1,res2)).sum()

        if res1 != res2:
            # freq_dict is a dict where values is a list: second element is the freq
            # need to multiply by 2 because need to count prob of (AB) AND (BA)
            n_exp_res1_res2_pairs = n_pairs_total * 0.05**2 *2
        elif res1 == res2:
            n_exp_res1_res2_pairs = n_pairs_total * 0.05**2

        scores[(res1, res2)] = -np.log10(n_obs_res1_res2_pairs / n_exp_res1_res2_pairs)
        obs_exp[(res1, res2)] = [n_obs_res1_res2_pairs, n_exp_res1_res2_pairs]
        print(res1, res2, obs_exp[(res1, res2)])

    return scores, obs_exp


def num_pairs(x):
    # return len(list(itertools.combinations(sorted(list(x)), 2)))
    # use permutations because
    # using pair[0] and pair[1] in lines 33,36
    return len(list(itertools.combinations(sorted(list(x)), 2)))

def to_pairs(x):
    return list(itertools.combinations(sorted(list(x)), 2))

def get_resnames(row):
    if row['resname_vdm'] in set(resname_dict.keys()):
        return True
    else:
        return False

######## START CODE HERE ###########
from sys import argv 
import sys
import pandas as pd, pickle as pkl
sys.path.append('/home/gpu/Sophia/combs/src/')
from combs import analysis

script, threepfive_csv_dir, ifg = argv
pdbinfo = pd.read_csv(threepfive_csv_dir+'%s_vdm_pdb_info.csv' %ifg)
ifginfo = pd.read_csv(threepfive_csv_dir+'%s_ifg_pdb_info.csv' %ifg)

merged = pd.merge(pdbinfo, ifginfo, on='iFG_count', suffixes=('_vdm', '_ifg'))
df = merged[['iFG_count', 'pdb', 'resname_ifg', 'resnum_ifg', 'chid_ifg', 'vdM_count', 'resname_vdm', 'resnum_vdm', 'chid_vdm', 'sequence_vdm', 'sec_struct_dssp_vdm']]
print(len(df), 'full df')
df = analysis.Analysis.remove_repeat_proteins(df)
print(len(df), 'removed repeats df')
pkl.dump(df, open('%s_3.5_norepeats.pkl'%ifg,'wb'))

#df = pkl.load(open('%s_3.5_norepeats.pkl'%ifg,'rb'))
corr = correlation(df)
pkl.dump(corr, open('noskip_%s_correlation.pkl'%ifg,'wb'))
