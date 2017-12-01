__all__ = ['make_freqaai_df']

from .. import cluster
import pickle as pkl, pandas as pd

def make_freqaai_df(skip10_pkl):   
    '''
    Inputs:
    Returns: dict of dictionaries. (1) counts for bb vdms. 
    (2) its freq. (3) counts for sc vdms. (4) its freq. 
    Note: some vdms will be counted as both bb and sc.'''
    
    bb_atoms = ['N', 'CA', 'C', 'O', 'OXT']
    df = pkl.load(open(skip10_pkl,'rb'))
    df = df[['resname_vdm', 'atom_names_vdm', 'dist_info']]

    # get bb and sc vdms that are within 3.5A
    bb = df[df.apply(cluster.Interactamer.has_bb_or_sc, bb_or_sc='bb', threepfive=True,axis=1)]
    sc = df[df.apply(cluster.Interactamer.has_bb_or_sc, bb_or_sc='sc', threepfive=True,axis=1)]
    num_within_threefive = len(set(list(pd.concat([bb, sc]).index.values)))
    print((num_within_threefive), 'num vdMs within 3.5A')
    print(len(bb), 'num bb interactamers within 3.5A')
    print(len(sc), 'num sc interactamers within 3.5A')
    bbcounts = pd.DataFrame(bb['resname_vdm'].value_counts())
    sccounts = pd.DataFrame(sc['resname_vdm'].value_counts())
    aai_df = pd.merge(sccounts,bbcounts, right_index=True, left_index=True, suffixes=('_sc', '_bb'))
    
    # drop rare AAs 
    for ix, row in aai_df.iterrows():
        rare = ['MSE', 'SEP', 'TPO', 'CSO']
        if ix in rare:
            aai_df.drop(ix, axis=0,inplace=True)
    
    # get frequencies from the counts
    for col in aai_df.columns.values:
        interaction_type = col.split('_')[2]
        total = sum(aai_df[col])
        freq_col = pd.Series(aai_df[col]/total, name='sdf')
        aai_df = pd.merge(aai_df, pd.DataFrame(freq_col), right_index=True, left_index=True)
        aai_df = aai_df.rename(index=str, columns={'sdf':'vdm_freq_'+interaction_type})


    return aai_df

