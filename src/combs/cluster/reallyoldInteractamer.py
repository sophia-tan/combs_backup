__all__ = ['make_all_rel_vdms', 'get_centroids', 'make_all_rel_vdms_hbond', 'make_all_rel_vdms_hbond_partial_ifg', 'unnamed']

import prody as pr
import pandas as pd
import numpy as np
import traceback
import os
import sys
import pickle
from ..apps.constants import interactamer_atoms, flip_sets, residue_sc_names
from ..analysis import Analyze
from ..analysis.cluster import cluster_adj_mat
from sklearn.neighbors import NearestNeighbors
import shutil

# class Interactamer():
#     def __init__(self):
#         pass

def listdir(path):
    return [file for file in os.listdir(path) if file[0] != '.']

def get_centroids(picklepath, radius=0.2):
    """This will cluster all vdms by iFG location (backbone rel_vdms) or by sidechain + iFG location (sc rel_vdms)
    and output new pickle files to a directory picklepath/clustered."""
    if picklepath[-1] != '/':
        picklepath += '/'
    if os.path.isdir(picklepath):
        for pickletype in listdir(picklepath):
            if pickletype == 'PHI_PSI':
                for phipsi_type in listdir(picklepath + pickletype):
                    for picklefile in listdir(picklepath + pickletype + '/' + phipsi_type + '/pickle/'):
                        with open(picklepath + pickletype + '/' + phipsi_type + '/pickle/'
                                  + picklefile, 'rb') as infile:
                            pick = pickle.load(infile)
                            if len(pick.shape) == 1:
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick, outfile)
                            else:
                                ifg_flat = [coords.flatten() for coords in pick[:, -2]]
                                nbrs = NearestNeighbors(metric='euclidean', radius=radius)
                                nbrs.fit(ifg_flat)
                                adj_mat = nbrs.radius_neighbors_graph(ifg_flat)
                                mems, cents = cluster_adj_mat(adj_mat)
                                outpath = picklepath + 'clustered/' + pickletype + '/' + phipsi_type + '/pickle/'
                                try:
                                    os.makedirs(outpath)
                                except:
                                    pass
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick[cents, :], outfile)
            else:
                for picklefile in listdir(picklepath + pickletype + '/pickle/'):
                    with open(picklepath + pickletype + '/pickle/' + picklefile, 'rb') as infile:
                        pick = pickle.load(infile)
                        outpath = picklepath + 'clustered/' + pickletype + '/pickle/'
                        try:
                            os.makedirs(outpath)
                        except:
                            pass
                        if len(pick.shape) == 1:
                            with open(outpath + picklefile, 'wb') as outfile:
                                pickle.dump(pick, outfile)
                        else:
                            if pickletype == 'SC':
                                sc_flat = [coords.flatten() for coords in pick[:, -3]]
                                ifg_flat = [coords.flatten() for coords in pick[:, -2]]
                                sc_ifg_flat = np.hstack((sc_flat, ifg_flat))
                                nbrs = NearestNeighbors(metric='euclidean', radius=radius)
                                nbrs.fit(sc_ifg_flat)
                                adj_mat = nbrs.radius_neighbors_graph(sc_ifg_flat)
                                mems, cents = cluster_adj_mat(adj_mat)
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick[cents, :], outfile)
                            else:
                                ifg_flat = [coords.flatten() for coords in pick[:, -2]]
                                nbrs = NearestNeighbors(metric='euclidean', radius=radius)
                                nbrs.fit(ifg_flat)
                                adj_mat = nbrs.radius_neighbors_graph(ifg_flat)
                                mems, cents = cluster_adj_mat(adj_mat)
                                with open(outpath + picklefile, 'wb') as outfile:
                                    pickle.dump(pick[cents, :], outfile)


def make_all_rel_vdms(path_to_csv, outpath, comb, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    dist_vdms_sc = get_vdms_sc(dist_vdms)
    make_rel_vdms(dist_vdms_sc, an, outpath + 'SC', comb, 'CA', 'N', 'C')
    dist_vdms_bb_N_CA = get_vdms_bb_N_CA(dist_vdms)
    make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    dist_vdms_bb_C_O = get_vdms_bb_C_O(dist_vdms)
    make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')
    dist_vdms_bb_phi_psi = get_vdms_bb_phi_psi(dist_vdms)
    for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
            outpath_phipsi = outpath + 'PHI_PSI/' + str(phi_low) + '_' + str(phi_high) + '_' + str(psi_low) \
                     + '_' + str(psi_high)
            def phi_psi_select(row):
                if type(row['sec_struct_phi_psi_vdm']) is float:
                    return False
                if 'None' in row['sec_struct_phi_psi_vdm']:
                    return False
                spltrow = row['sec_struct_phi_psi_vdm'].split()
                if (phi_low <= float(spltrow[0]) < phi_high) and (psi_low <= float(spltrow[1]) < psi_high):
                    return True
                else:
                    return False
            dist_vdms_phipsi = dist_vdms_bb_phi_psi[dist_vdms_bb_phi_psi.apply(phi_psi_select, axis=1)]
            make_rel_vdms(dist_vdms_phipsi, an, outpath_phipsi, comb, 'CA', 'N', 'C')


def get_sc_hb(row):
    try:
        index = [i for i, rrn in enumerate([n for n in row['rel_resnums'] if n != '-']) if rrn == '0']
        return any([any(
            {y for names in [row['vdM_atom_names'].strip('()').split(') (')[ind]] for y in names.split() if
             y not in ['N', 'O', 'H', 'OXT', 'H1', 'H2', 'H3']}) for ind in index])
    except:
        return False


def get_bb_hb(row):
    try:
        index = [i for i, rrn in enumerate([n for n in row['rel_resnums'] if n != '-']) if rrn == '0']
        return any([any(
            {y for names in [row['vdM_atom_names'].strip('()').split(') (')[ind]] for y in names.split() if
             y in ['N', 'O', 'H']}) for ind in index])
    except:
        return False


def make_all_rel_vdms_hbond(path_to_csv, outpath, comb, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm['rel_resnums'].str.contains('0')]
    hbond_vdms_sc = hbond_vdms[hbond_vdms.apply(get_sc_hb, axis=1)]
    hbond_vdms_bb = hbond_vdms[hbond_vdms.apply(get_bb_hb, axis=1)]
    df_dist_vdms_sc = pd.merge(dist_vdms, hbond_vdms_sc, on=['iFG_count', 'vdM_count'])
    df_dist_vdms_bb = pd.merge(dist_vdms, hbond_vdms_bb, on=['iFG_count', 'vdM_count'])
    dist_vdms_sc = get_vdms_sc(df_dist_vdms_sc)
    make_rel_vdms(dist_vdms_sc, an, outpath + 'SC', comb, 'CA', 'N', 'C')
    dist_vdms_bb_N_CA = get_vdms_bb_N_CA(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    dist_vdms_bb_C_O = get_vdms_bb_C_O(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')
    dist_vdms_bb_phi_psi = get_vdms_bb_phi_psi(df_dist_vdms_bb)
    for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
            outpath_phipsi = outpath + 'PHI_PSI/' + str(phi_low) + '_' + str(phi_high) + '_' + str(psi_low) \
                     + '_' + str(psi_high)
            def phi_psi_select(row):
                if type(row['sec_struct_phi_psi_vdm']) is float:
                    return False
                if 'None' in row['sec_struct_phi_psi_vdm']:
                    return False
                spltrow = row['sec_struct_phi_psi_vdm'].split()
                if (phi_low <= float(spltrow[0]) < phi_high) and (psi_low <= float(spltrow[1]) < psi_high):
                    return True
                else:
                    return False
            dist_vdms_phipsi = dist_vdms_bb_phi_psi[dist_vdms_bb_phi_psi.apply(phi_psi_select, axis=1)]
            make_rel_vdms(dist_vdms_phipsi, an, outpath_phipsi, comb, 'CA', 'N', 'C')


def make_all_rel_vdms_hbond_partial_ifg(path_to_csv, outpath, comb, ifg_name=None, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    if ifg_name:
        def func(row):
            if (ifg_name in row['iFG_atom_names']) and ('0' in row['rel_resnums']):
                return True
            else:
                return False
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm.apply(func, axis=1)]
    else:
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm['rel_resnums'].str.contains('0')]
    hbond_vdms_sc = hbond_vdms[hbond_vdms.apply(get_sc_hb, axis=1)]
    hbond_vdms_bb = hbond_vdms[hbond_vdms.apply(get_bb_hb, axis=1)]
    df_dist_vdms_sc = pd.merge(dist_vdms, hbond_vdms_sc, on=['iFG_count', 'vdM_count'])
    df_dist_vdms_bb = pd.merge(dist_vdms, hbond_vdms_bb, on=['iFG_count', 'vdM_count'])
    dist_vdms_sc = get_vdms_sc(df_dist_vdms_sc)
    make_rel_vdms(dist_vdms_sc, an, outpath + 'SC', comb, 'CA', 'N', 'C')
    dist_vdms_bb_N_CA = get_vdms_bb_N_CA(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    dist_vdms_bb_C_O = get_vdms_bb_C_O(df_dist_vdms_bb)
    make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')
    dist_vdms_bb_phi_psi = get_vdms_bb_phi_psi(df_dist_vdms_bb)
    for phi_low, phi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
        for psi_low, psi_high in zip(range(-180, 180, 36), range(-144, 216, 36)):
            outpath_phipsi = outpath + 'PHI_PSI/' + str(phi_low) + '_' + str(phi_high) + '_' + str(psi_low) \
                     + '_' + str(psi_high)
            def phi_psi_select(row):
                if type(row['sec_struct_phi_psi_vdm']) is float:
                    return False
                if 'None' in row['sec_struct_phi_psi_vdm']:
                    return False
                spltrow = row['sec_struct_phi_psi_vdm'].split()
                if (phi_low <= float(spltrow[0]) < phi_high) and (psi_low <= float(spltrow[1]) < psi_high):
                    return True
                else:
                    return False
            dist_vdms_phipsi = dist_vdms_bb_phi_psi[dist_vdms_bb_phi_psi.apply(phi_psi_select, axis=1)]
            make_rel_vdms(dist_vdms_phipsi, an, outpath_phipsi, comb, 'CA', 'N', 'C')


def make_all_interactamers_hbond_partial_ifg(path_to_csv, outpath, comb, ifg_name=None, seq_distance=5):
    if outpath[-1] != '/':
        outpath += '/'
    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=seq_distance)
    if ifg_name:
        def func(row):
            if (ifg_name in row['iFG_atom_names']) and ('0' in row['rel_resnums']):
                return True
            else:
                return False
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm.apply(func, axis=1)]
    else:
        hbond_vdms = an.ifg_hbond_vdm[an.ifg_hbond_vdm['rel_resnums'].str.contains('0')]
    hbond_vdms_sc = hbond_vdms[hbond_vdms.apply(get_sc_hb, axis=1)]
    hbond_vdms_bb = hbond_vdms[hbond_vdms.apply(get_bb_hb, axis=1)]
    df_dist_vdms_sc = pd.merge(dist_vdms, hbond_vdms_sc, on=['iFG_count', 'vdM_count'])
    # df_dist_vdms_bb = pd.merge(dist_vdms, hbond_vdms_bb, on=['iFG_count', 'vdM_count'])
    dist_vdms_sc = get_vdms_sc(df_dist_vdms_sc)
    make_interactamers(dist_vdms_sc, an, outpath + 'SC', comb)




def make_rel_vdms(df, an, outpath, comb, origin_atom, plane_atom1, plane_atom2):
    if outpath[-1] != '/':
        outpath += '/'
    resns = set(df.groupby('resname_vdm').groups).intersection(set(interactamer_atoms.keys()))
    for resn in resns:
        vdms = parse_interactamers_aa(df, an, resn)
        if vdms:
            pdb_path = outpath + 'pdbs/' + resn + '/'
            try:
                os.makedirs(pdb_path)
            except:
                pass
            picklepath = outpath + 'pickle/'
            try:
                os.makedirs(picklepath)
            except:
                pass
            rel_vdm_output = []

            print('Making relative vdMs for ' + resn + '...')
            for vdm in vdms:
                try:
                    ifg_count = int(repr(vdm).split('_')[1])
                    vdm_count = int(repr(vdm).split('_')[3])
                    bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs = make_rel_vdm_coords(vdm, comb, origin_atom,
                                                                                  plane_atom1, plane_atom2)
                    for i, bbc, scc, sccon, scccs, ic, icon, iccs, pdb in zip(range(len(bb_coords)), bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs):
                        rel_vdm_output.append([ifg_count, vdm_count, i + 1, bbc, scc, sccon, scccs, ic, icon, iccs, resn])
                        string = repr(pdb).split()[1].split('_')
                        pr.writePDB(pdb_path + '_'.join(string[:-1]) + '_iFlip_' + str(i + 1)
                                    + '_' + string[-1] + '_oriented.pdb.gz', pdb)
                except Exception:
                    traceback.print_exc()
            # output format = [ifg_count, vdm_count, ifg_flip, bb_coords, sc_coords_all, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, vdm_resn]

            if rel_vdm_output:
                with open(picklepath + resn + '_rel_vdms.pickle', 'wb') as f:
                    pickle.dump(np.array(rel_vdm_output, dtype=object), f)
            else:
                shutil.rmtree(pdb_path)


def make_interactamers(df, an, outpath, comb):
    if outpath[-1] != '/':
        outpath += '/'
    resns = set(df.groupby('resname_vdm').groups).intersection(set(interactamer_atoms.keys()))
    for resn in resns:
        print(resn)
        for key in interactamer_atoms[resn].keys():
            origin_atom, plane_atom1, plane_atom2 = interactamer_atoms[resn][key]
            vdms = parse_interactamers_aa(df, an, resn)
            if vdms:
                pdb_path = outpath + 'vdM/' + resn + '/' + key +'/'
                try:
                    os.makedirs(pdb_path)
                except:
                    pass
                picklepath = outpath + 'pickle/'
                try:
                    os.makedirs(picklepath)
                except:
                    pass
                rel_vdm_output = []

                print('Making interactamer vdMs for ' + resn + ', ' + key + '...')
                for vdm in vdms:
                    try:
                        ifg_count = int(repr(vdm).split('_')[1])
                        vdm_count = int(repr(vdm).split('_')[3])
                        bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs = make_rel_vdm_coords(vdm, comb, origin_atom,
                                                                                      plane_atom1, plane_atom2)
                        for i, bbc, scc, sccon, scccs, ic, icon, iccs, pdb in zip(range(len(bb_coords)), bb_coords, sc_coords, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbs):
                            rel_vdm_output.append([ifg_count, vdm_count, i + 1, bbc, scc, sccon, scccs, ic, icon, iccs, resn])
                            string = repr(pdb).split()[1].split('_')
                            pr.writePDB(pdb_path + '_'.join(string[:-1]) + '_iFlip_' + str(i + 1)
                                        + '_' + string[-1] + '_oriented.pdb.gz', pdb)
                    except Exception:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        print(exc_type, fname, exc_tb.tb_lineno, exc_obj)
                # output format = [ifg_count, vdm_count, ifg_flip, bb_coords, sc_coords_all, sc_coords_ON, sc_coords_CS, ifg_coords, ifg_ON_coords, ifg_CS_coords, vdm_resn]

                if rel_vdm_output:
                    with open(picklepath + resn + '_rel_vdms.pickle', 'wb') as f:
                        pickle.dump(np.array(rel_vdm_output, dtype=object), f)
                else:
                    shutil.rmtree(pdb_path)


def get_vdms_sc(df):
    def sc(row):
        if set(row['atom_names_vdm'].split()) - {'N', 'C', 'CA', 'O', 'OXT'} != set():
            return True
        else:
            return False
    return df[df.apply(sc, axis=1)]


def get_vdms_bb_N_CA(df):
    def bb(row):
        if any(name in row['atom_names_vdm'].split() for name in ['N', 'CA']):
            if all(name not in row['atom_names_vdm'].split() for name in ['C', 'O']):
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]


def get_vdms_bb(df, NCAorCO):
    # added by Sophia. NCAorCO can be 'NCA' or 'CO'
    def bb(row):
        if NCAorCO=='NCA':
            if any(name in row['atom_names_vdm'].split() for name in ['N', 'CA']):
                return True
            else:
                return False
        elif NCAorCO=='CO':
            if any(name in row['atom_names_vdm'].split() for name in ['C', 'O']):
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]


def get_vdms_bb_C_O(df):
    def bb(row):
        if any(name in row['atom_names_vdm'].split() for name in ['C', 'O']):
            if 'N' not in row['atom_names_vdm'].split():
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]


def get_vdms_bb_phi_psi(df):
    def bb(row):
        if 'N' in row['atom_names_vdm'].split():
            if any(name in row['atom_names_vdm'].split() for name in ['C', 'O']):
                return True
            else:
                return False
        else:
            return False
    return df[df.apply(bb, axis=1)]


def parse_interactamers_aa(df, an, resn):
    #use dist_vdms as dataframe (df)
    df = df[df['resname_vdm'] == resn]
    return an.parse_vdms(df)


def make_rel_vdm_coords(pdb, comb, origin_atom, plane_atom1, plane_atom2):
    try:
        origin_coords = pdb.select('chain X and resnum 10 and name ' + origin_atom).getCoords()[0]
        pdb_coords = pdb.getCoords()
        pdb_coords_neworigin = pdb_coords - origin_coords
        plane_atom1_coords = pdb.select('chain X and resnum 10 and name ' + plane_atom1).getCoords()[0] - origin_coords
        plane_atom2_coords = pdb.select('chain X and resnum 10 and name ' + plane_atom2).getCoords()[0] - origin_coords
        x_norm = plane_atom1_coords / np.linalg.norm(plane_atom1_coords)
        orthvec = np.cross(plane_atom1_coords, plane_atom2_coords)
        z_norm = -1 * orthvec / np.linalg.norm(orthvec)
        orthvec2 = np.cross(plane_atom1_coords, orthvec)
        y_norm = orthvec2 / np.linalg.norm(orthvec2)
        R = np.array([x_norm, y_norm, z_norm])
        pdb_coords_neworigin_rot = np.dot(pdb_coords_neworigin, R.T)
        pdbcopy = pdb.copy()
        pdbcopy.setCoords(pdb_coords_neworigin_rot)
        new_origin_coords = pdbcopy.select('chain X and resnum 10 and name ' + origin_atom).getCoords()[0]
        new_plane_atom1_coords = pdbcopy.select('chain X and resnum 10 and name ' + plane_atom1).getCoords()[0]
        new_plane_atom2_coords = pdbcopy.select('chain X and resnum 10 and name ' + plane_atom2).getCoords()[0]
        # vdm_coords = [y for x in [new_origin_coords, new_plane_atom1_coords, new_plane_atom2_coords] for y in x]
        bb_coords = np.array([new_origin_coords, new_plane_atom1_coords, new_plane_atom2_coords])
        if pdbcopy.select('chain X and resnum 10 and name CA').getResnames()[0] != 'GLY':
            # sc_coords = pdbcopy.select('chain X and resnum 10 and sidechain and not element H D').getCoords()
            sc_coords = make_sc_atom_coords(pdbcopy)
            sc_coords_ON = make_sc_atom_coords_ON(pdbcopy)
            sc_coords_CS = make_sc_atom_coords_CS(pdbcopy)
        else:
            sc_coords = None
            sc_coords_ON = None
            sc_coords_CS = None
        ifg_sels = make_ifg_atom_sele(pdbcopy, comb)
        ifg_coords = []
        ifg_CS_coords = []
        ifg_ON_coords = []
        pdbcopies = []
        orig_sel = ifg_sels[0]
        for ifg_sel in ifg_sels:
            # ifg_coords.append([coord for ifg in ifg_sel for sel in ifg for coord in sel.getCoords()])
            coords = np.array([sel.getCoords() for ifg in ifg_sel for sel in ifg])
            coords_ON = [ifg.select('element O N').getCoords() for ifg in ifg_sel
                                  if ifg.select('element O N') is not None][0]
            coords_CS = [ifg.select('element C S').getCoords() for ifg in ifg_sel
                                  if ifg.select('element C S') is not None][0]
            ifg_coords.append(coords)
            ifg_CS_coords.append(coords_CS)
            ifg_ON_coords.append(coords_ON)
            for sel, coor in zip(orig_sel, coords):
                sel.setCoords(coor)
            pdbcopies.append(pdbcopy.copy())
        return [bb_coords]*len(ifg_coords), [sc_coords]*len(ifg_coords), [sc_coords_ON]*len(ifg_coords), \
               [sc_coords_CS]*len(ifg_coords), ifg_coords, ifg_ON_coords, ifg_CS_coords, pdbcopies
    except Exception:
        traceback.print_exc()
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno, exc_obj)



def make_sc_atom_coords(pdb):
    resname = pdb.select('chain X and resnum 10 and name CA').getResnames()[0]
    sc_coords = np.array([pdb.select('chain X and resnum 10 and name ' + name).getCoords()[0]
                          for name in residue_sc_names[resname]])
    return sc_coords

def make_sc_atom_coords_ON(pdb):
    resname = pdb.select('chain X and resnum 10 and name CA').getResnames()[0]
    sc_coords = np.array([pdb.select('chain X and resnum 10 and name ' + name).getCoords()[0]
                          for name in residue_sc_names[resname] if name[0] in {'N', 'O'}])
    return sc_coords

def make_sc_atom_coords_CS(pdb):
    resname = pdb.select('chain X and resnum 10 and name CA').getResnames()[0]
    sc_coords = np.array([pdb.select('chain X and resnum 10 and name ' + name).getCoords()[0]
                          for name in residue_sc_names[resname] if name[0] in {'C', 'S'}])
    return sc_coords

# This code is problematic because what if the iFG has a PHE or TYR in it, where more than 2 atoms need to be flipped at once?
def make_ifg_atom_sele(pdb, comb):
    """uses iFG definitions in comb object to select iFGs in the parsed protein object that have all atoms
    and occupancies = 1.
    """
    # There is a problem with this code: What if one wants to select atoms from a HEME as an iFG?
    # It is not represented in the one_letter_code dictionary...
    ifgs = []
    if comb.num_res_ifg == 1:
        poss_ifg_sel = pdb.select('chain Y and resnum 10')
        if poss_ifg_sel is not None:
            ifg_resindices, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames = poss_ifg_sel.getResnames()[indices]
            for ifg_resindex, ifg_resname in zip(ifg_resindices, ifg_resnames):
                    ifg_selection = pdb.select('resindex ' + str(ifg_resindex) + ' and name '
                                                          + comb.ifg_sele_dict[1][ifg_resname])
                    if ifg_selection is not None:
                        num_atoms = len(ifg_selection)
                        name_list = comb.ifg_sele_dict[1][ifg_resname].split()
                        name_set = set(name_list)
                        if num_atoms == len(name_set):
                            if all(ifg_selection.getResnums() > 0):
                                ifgs.append([ifg_selection.select('name ' + name) for name in name_list])
                            for flip_set in flip_sets:
                                if flip_set.issubset(name_set):
                                    new_names = np.array(name_list)
                                    flip_names = list(flip_set)
                                    i, = np.where(new_names == flip_names[0])
                                    j, = np.where(new_names == flip_names[1])
                                    temp = new_names[i]
                                    new_names[i] = new_names[j]
                                    new_names[j] = temp
                                    ifgs.append([ifg_selection.select('name ' + name) for name in new_names])
        return ifgs
    else:
        poss_ifg_sel = pdb.select('chain Y and resnum 10to11')
        if poss_ifg_sel is not None:
            ifg_resindices_cat_list, indices = np.unique(poss_ifg_sel.getResindices(), return_index=True)
            ifg_resnames_cat_list = poss_ifg_sel.getResnames()[indices]
            ifg_resindex_pairs = [ifg_resindices_cat_list[i:i + 2] for i in range(0, len(ifg_resindices_cat_list), 2)]
            ifg_resname_pairs = [ifg_resnames_cat_list[i:i + 2] for i in range(0, len(ifg_resnames_cat_list), 2)]
            for ifg_resindex_pair, ifg_resname_pair in zip(ifg_resindex_pairs, ifg_resname_pairs):
                resind1, resind2 = ifg_resindex_pair
                resname1, resname2 = ifg_resname_pair
                try:
                    ifg_selection1 = pdb.select('(resindex ' + str(resind1) + ' and name '
                                                      + comb.ifg_sele_dict[1][resname1]+')'
                                                      )
                    ifg_selection2 = pdb.select('(resindex ' + str(resind2) + ' and name '
                                                + comb.ifg_sele_dict[2][resname2] + ')')
                except KeyError:
                    print('Non-canonical residue in iFG, skipping.')
                    ifg_selection1 = None
                if ifg_selection1 is not None and ifg_selection2 is not None:
                    num_atoms = len(ifg_selection1 | ifg_selection2)
                    name_list1 = comb.ifg_sele_dict[1][resname1].split()
                    name_list2 = comb.ifg_sele_dict[2][resname2].split()
                    name_set1 = set(name_list1)
                    name_set2 = set(name_list2)
                    if num_atoms == (len(name_set1) + len(name_set2)):
                        if all(ifg_selection1.getResnums() > 0) and all(ifg_selection2.getResnums() > 0):
                            temp_ifg = [ifg_selection1.select('name ' + name) for name in name_list1]
                            temp_ifg.extend(ifg_selection2.select('name ' + name) for name in name_list2)
                            ifgs.append(temp_ifg)
                            flip1 = False
                            flip2 = False
                            for flip_set in flip_sets:
                                if flip_set.issubset(name_list1):
                                    flip1 = True
                                    new_names1 = np.array(name_list1)
                                    flip_names = list(flip_set)
                                    i, = np.where(new_names1 == flip_names[0])
                                    j, = np.where(new_names1 == flip_names[1])
                                    temp = new_names1[i]
                                    new_names1[i] = new_names1[j]
                                    new_names1[j] = temp
                                    temp_ifg = [ifg_selection1.select('name ' + name) for name in new_names1]
                                    temp_ifg.extend(ifg_selection2.select('name ' + name) for name in name_list2)
                                    ifgs.append(temp_ifg)
                                if flip_set.issubset(name_list2):
                                    flip2 = True
                                    new_names2 = np.array(name_list2)
                                    flip_names = list(flip_set)
                                    i, = np.where(new_names2 == flip_names[0])
                                    j, = np.where(new_names2 == flip_names[1])
                                    temp = new_names1[i]
                                    new_names1[i] = new_names1[j]
                                    new_names1[j] = temp
                                    temp_ifg = [ifg_selection1.select('name ' + name) for name in name_list1]
                                    temp_ifg.extend(ifg_selection2.select('name ' + name) for name in new_names2)
                                    ifgs.append(temp_ifg)
                            if flip1 and flip2:
                                temp_ifg = [ifg_selection1.select('name ' + name) for name in new_names1]
                                temp_ifg.extend(ifg_selection2.select('name ' + name) for name in new_names2)
                                ifgs.append(temp_ifg)
        return ifgs

### SOPHIA ADDED ###
def has_bb_or_sc(row, bb_or_sc,threepfive=False):
    # determine whether there are bb or sc vdmatoms
    # threepfive is an option for calculating freq_aai where you
    # only want to get the bb & sc vdMs if the interacting atoms
    # are within 3.5A
    bb_atoms = ['N', 'CA', 'C', 'O', 'OXT']
    vdmatoms = row['atom_names_vdm'].split(' ')
    num_bb_in_vdm = sum([x in bb_atoms for x in vdmatoms])
    has_bb = num_bb_in_vdm > 0
    has_sc = len(vdmatoms) > num_bb_in_vdm
    bb_bool = np.array([x in bb_atoms for x in vdmatoms]) 
    if bb_or_sc=='bb':
        if threepfive==True:
            close = within_threepfive(row, bb_atoms,bbsc='bb')
            return has_bb and close
        return has_bb 
    elif bb_or_sc=='sc':
        if threepfive==True:
            close = within_threepfive(row, bb_atoms,bbsc='sc')
            return has_sc and close
        return has_sc 

def within_threepfive(row,bb_atoms,bbsc):
    # returns whether or not the interacting atom is within 3.5A
    dist = row['dist_info']
    dist = dist.strip('()').split(') (')
    dist = [x.split(' ') for x in dist] 
    # this is now a list of lists where the inner list is 
    # element0: iFG, element1: vdM, element2: distance
    # check if the vdMs are bb or sc atoms
    if bbsc=='bb':
        vdms = [v for v in dist if v[1] in bb_atoms]
    elif bbsc=='sc':
        vdms = [v for v in dist if v[1] not in bb_atoms]
    distances = [float(v[2]) for v in vdms]
    numclose = sum([d<=3.5 for d in distances])
    return numclose > 0

def unnamed(ifg_dir, comb):
    path_to_csv = ifg_dir+'/csv/'
    outpath = ifg_dir+'/clusters/'

    an = Analyze(path_to_csv)
    dist_vdms = an.get_distant_vdms(seq_distance=10)
    vdms_bb = dist_vdms[dist_vdms.apply(has_bb_or_sc, bb_or_sc='bb',axis=1)]
    vdms_sc = dist_vdms[dist_vdms.apply(has_bb_or_sc, bb_or_sc='sc',axis=1)]
    make_interactamers(vdms_sc, an, outpath + 'transformed', comb)

    dist_vdms_bb_N_CA = get_vdms_bb(dist_vdms,NCAorCO='NCA')
    dist_vdms_bb_C_O = get_vdms_bb(dist_vdms,NCAorCO='CO')
    #make_rel_vdms(dist_vdms_bb_N_CA, an, outpath + 'N_CA', comb, 'N', 'H', 'CA')
    #make_rel_vdms(dist_vdms_bb_C_O, an, outpath + 'C_O', comb, 'C', 'O', 'CA')
