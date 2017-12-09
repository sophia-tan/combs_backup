from prody import *
import pandas as pd

def append_phi_psi(phipsi, parsed_res):
    try:
        phi = calcPhi(parsed_res)
    except Exception:
        phi = None
    try:
        psi = calcPsi(parsed_res)
    except:
        psi = None
    if phi is not None and psi is not None:
        phipsi.append([phi, psi])
    elif phi is None and psi is not None:
        phipsi.append([None, psi])
    elif phi is not None and psi is None:
        phipsi.append([phi, None])
    return phipsi

def get_phi_psi_list(pdb, chain):
        phipsi = []

        parsed = parsePDB(pdb)
        ca = parsed.select('chain %s and calpha' % chain) # get ca for resnums
        resnums = ca.getResnums()
        for resn in resnums:
            parsed_res = parsed[chain, resn]
            phipsi = append_phi_psi(phipsi, parsed_res)
        return phipsi

def get_ss_string(frag):
    phipsi_dict = {}
    phipsi_dict['H'] = [[-180, 0], [-100,45]] # [phi, psi]
    phipsi_dict['S'] = [[-180,-45], [45,225]]
    phipsi_dict['T'] = [[0, 180], [-90,90]] # this is actually called "turn" in the Hovmoller 2002 paper

    ss_list = '' 
    for ix, res in enumerate(frag):
        assigned = 0
        phi, psi = res[0], res[1]
        for k, v in phipsi_dict.items():
            try:
                if phi >= v[0][0] and phi <= v[0][1]:
                    if psi >= v[1][0] and psi <= v[1][1]:
                        ss_list += k
                        assigned += 1

            except: # will catch None
                if phi == None:
                  if psi >= v[1][0] and psi <= v[1][1]:
                      ss_list += k
                      assigned += 1
                if psi == None:
                    if phi >= v[0][0] and phi <= v[0][1]:
                        ss_list += k
                        assigned += 1
                else:
                    pass
        # if ss wasn't added to ss_list, mark as '-'
        if assigned == 0:
            ss_list += '-'
    return ss_list 

def generate_ss_df():

    col_names = ['all helix', 'Nterm strand-helix', 'helix-Cterm turn', 'helix-Cterm strand', 'all strand', 'all loop', 'no ss', '% no ss']
    
    ss_df = pd.DataFrame(index=[1,2,3,4], columns=col_names)
    
    # fill with zeros
    for ix, row in ss_df.iterrows():
        for col in col_names:
            row[col] = 0
    return ss_df


def inc_ss_df(df, frag, num):
    ss_string = frag

    if abs(num) == 1:
        if ss_string == 'HH':
            ss = 'all helix'
        if ss_string == 'SH':
            ss = 'Nterm strand-helix'
        if ss_string == 'HT':
            ss = 'helix-Cterm turn'
        if ss_string == 'HS':
            ss = 'helix-Cterm strand'
        if ss_string == 'SS':
            ss = 'all strand'
        if ss_string == 'TT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            df.ix[num, 'no ss'] += 1
        
    elif abs(num) == 2:
        if ss_string == 'HHH':
            ss = 'all helix'
        if ss_string == 'SHH':
            ss = 'Nterm strand-helix'
        if ss_string == 'HHT':
            ss = 'helix-Cterm turn'
        if ss_string == 'HHS':
            ss = 'helix-Cterm strand'
        if ss_string == 'SSS':
            ss = 'all strand'
        if ss_string == 'TTT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            df.ix[num, 'no ss'] += 1

    elif abs(num) == 3:
        if ss_string == 'HHHH':
            ss = 'all helix'
        if ss_string == 'SHHH':
            ss = 'Nterm strand-helix'
        if ss_string == 'HHHT':
            ss = 'helix-Cterm turn'
        if ss_string == 'HHHS':
            ss = 'helix-Cterm strand'
        if ss_string == 'SSSS':
            ss = 'all strand'
        if ss_string == 'TTTT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            df.ix[num, 'no ss'] += 1
            
    elif abs(num) == 4 or abs(num) == 5: # treat 4 and 5 the same
        if ss_string[:5] == 'HHHHH':
            ss = 'all helix'
        if ss_string[:5] == 'SHHHH':
            ss = 'Nterm strand-helix'
        if ss_string[:5] == 'HHHHT':
            ss = 'helix-Cterm turn'
        if ss_string[:5] == 'HHHHS':
            ss = 'helix-Cterm strand'
        if ss_string[:5] == 'SSSSS':
            ss = 'all strand'
        if ss_string[:5] == 'TTTTT':
            ss = 'all loop'
        try:
            df.ix[num, ss] += 1
        except:
            df.ix[num, 'no ss'] += 1
        
    else:
        print(ss_string)
    return df

def perc_no_ss(ss_df):
    # count how many aren't accounted for in terms of ss
    for ix, row in ss_df.iterrows():
        total = sum(row)
    
        ss_df.ix[ix, '% no ss'] = row['no ss']/total* 100
    return ss_df 
