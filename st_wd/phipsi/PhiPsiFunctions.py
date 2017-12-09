import numpy as np
import pandas as pd
import math
import pickle as pkl

def rounddown(x): 
    # round down to nearest 10
    return int(math.floor(x/10))*10

def create_phi_psi_matrix():
    # in matrix, phi is x axis, psi is y axis!
    ls = np.arange(-180, 190, 10) # these numbers are the LOWER BOUNDS. ex) for the -180 bin, numbers that are -180 <= x < -180+10 is fine
    df = pd.DataFrame(index=ls, columns=ls)
    for x in ls:
        for y in ls:
            df.ix[x,y] = 0
    return df

def inc_phi_psi_df(counts_df, vdm_df):
    phi_psi = vdm_df['sec_struct_phi_psi_vdm']
    for row in phi_psi:
        if type(row) != float: # if type(row) == float, it is NECESSARILY (i checked) 'np.nan'
            row = row.split(' ') # checked that len(row) == 2
            if 'None' not in row: # n-term or c-term residue, makes up 1% of vdms
                row = [float(x) for x in row]
                # phi = columns, psi = index
                phi = rounddown(row[0])
                psi = rounddown(row[1])
                try:
                    counts_df.ix[psi, phi] += 1 # reversed phi/psi bc col/row!
                except:
                    print(row, phi, psi)
    return counts_df

def inc_phi_psi_df_for_db(counts_df, phi, psi):
    if phi != None and psi != None:
        # phi = columns, psi = index
        phi = rounddown(phi)
        psi = rounddown(psi)
        try:
            counts_df.ix[psi, phi] += 1 # reversed phi/psi bc col/row!
        except:
            print('couldn\'t add to df', phi, psi)
    return counts_df

def log(x):
    return pd.Series([np.log(i) for i in x], name=x.name)

def normalize(df):
    # calculate total counts in df
    sums = (sum(df.values))
    total = sum(sums)
    for ix, row in df.iterrows():
        for col in df.columns.values:
            df.ix[ix, col] = df.ix[ix,col] / total*100
    return df


def log_norm_df(df_path):

    df = pd.read_pickle(df_path)
    df = df[df.columns].astype(float)
    df = df.sort_index(ascending=False)
    
    # turn 0's to 1's for log transformation
    for ix, row in df.iterrows():
        for col in df.columns.values:
            if df.ix[ix, col] == 0:
                df.ix[ix,col] = 1
    
    df = df.apply(log)
    # normalize whole plot so that each bin represents the %
    df = normalize(df)
    
    return df
