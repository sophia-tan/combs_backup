__all__ = ['hierarchical_cluster', 'do_hier_clus_for_bb_and_sc']

from prody import *
import numpy as np
import sys
import os, pickle as pkl
import shutil
from ..apps.constants import interactamer_atoms
import scipy
from scipy.cluster import hierarchy as hier
from scipy.cluster.hierarchy import fcluster
import itertools
from ..analysis.cluster import cluster_adj_mat
from sklearn.neighbors import NearestNeighbors, KNeighborsClassifier

def concat_hier_pca_inputs(list_iFGcoords, list_pdbs):
    '''Concatenates inputs for hierarchical clustering and PCA. Use if you're analyzing >1 VDM AA 
    (ex: ASN+GLN or ASP+GLU)'''
    multiple_iFGcoords = np.vstack(list_iFGcoords)
    if len(list_pdbs)==1:
        multiple_pdbs = list_pdbs[0]
    elif len(list_pdbs)==2:
        multiple_pdbs=list_pdbs[0]+list_pdbs[1]
    return multiple_iFGcoords, multiple_pdbs

def getpdbpath(pkldata,ifg,pkldir,BBorSC):
    '''helper function. gets pdbpath from pickle file'''
    ls = []
    for item in pkldata:
        ifgcount,vdmcount,iflip,vflip = item[0],item[1],item[2],item[3]
        resn = item[11]
        resngroup = list(interactamer_atoms[resn].keys())[0]
        path = 'iFG_%s_vdM_%s_iFlip_%s_vFlip_%s_%s_oriented.pdb.gz'\
            % (ifgcount,vdmcount,iflip,vflip,ifg)
        if BBorSC == 'SC':
            pdbdir = pkldir+'pdbs/%s/%s/'%(resn,resngroup)
        elif BBorSC == 'N_CA' or BBorSC=='C_O':
            pdbdir = pkldir+'pdbs/%s/'%(resn) # doesn't have resngroup folder
        path = pdbdir+path
        ls.append([[resn, ifgcount, vdmcount], path])
    return ls

def do_hier_clus_for_bb_and_sc(combsoutputdir, ifg, vdmresidues, print_clus):
    '''print_clus needs to be boolean'''
    # iterate for bb and sc interactions
    BBorSC = ['SC']
    #BBorSC = ['SC', 'N_CA', 'C_O'] 
    vdmresidues = vdmresidues.split(',')
    vdmnames = '+'.join(vdmresidues) # for naming pickle file
    for typ in BBorSC:
        # make outputdir 
        outdir = combsoutputdir+ifg+'/clusters/%s/clustered/'%typ
        try:
            os.mkdir(outdir)
        except:
            pass
        # load ifgcoords and ifgcount/vdmcount info
        pkldir = combsoutputdir+ifg+'/clusters/%s/'%typ
        pkls = [pkl.load(open(pkldir+'pickle/%s_rel_vdms.pickle'%vdm,'rb')) for vdm in vdmresidues]
        ifgcoords = [x[:,8] for x in pkls]
        ifgcoords = [np.array([y.flatten() for y in x]) for x in ifgcoords]
        pdbpaths = [getpdbpath(x,ifg,pkldir,typ) for x in pkls]
        ifgcoords, pdbs = concat_hier_pca_inputs(ifgcoords,pdbpaths)

        if print_clus==False:
            sorted_clusters=hierarchical_cluster(ifgcoords, pdbs, max_d=0.7, print_clusters = False, outputdir = None)
        elif print_clus==True:
            sorted_clusters=hierarchical_cluster(ifgcoords, pdbs, max_d=0.7, print_clusters = True, outputdir = outdir,vdmnames=vdmnames,BBorSC=typ)
                
        pklfile = open(outdir+'clusterinfo_%s_%s_%s.pkl'%(ifg,vdmnames,typ), 'wb')
        pkl.dump(sorted_clusters, pklfile)

def make_rmsd_matrix(iFGcoords_array):
    pair_dist = scipy.spatial.distance.pdist(iFGcoords_array)
    pair_dist_matrix = scipy.spatial.distance.squareform(pair_dist)
    # pair_dist_matrix is the rmsd matrix
    return pair_dist_matrix
    
def hierarchical_cluster(iFGcoords, pdbs_of_iFG_coords, max_d=0.9, print_clusters = False, outputdir = None, vdmnames=None,BBorSC=None):
    '''Performs hierarchical clustering on a set of iFGcoords array
    Inputs: max_d (hierarchical clustering cutoff), ifgcoordsarray, pdbs. Option to print the 
            clusters to a directory to view what the clusters look like.
    Returns: sorted_clusters (list where each element is a list of all the pdbs in that cluster)
    '''

    rmsd_mat = make_rmsd_matrix(iFGcoords)
    '''
    Z = hier.linkage(iFGcoords,method='single')
    #Z = hier.linkage(scipy.spatial.distance.squareform(rmsd_mat),method='single')
    clusters = fcluster(Z, max_d, criterion='distance') # retrieve clusters
    #print(clusters)
    n_clus = max(clusters)
    unsorted_clusters = []
    for n in range(1,n_clus+1):
        members = [np.where(clusters == n)][0][0]
        members = [pdbs_of_iFG_coords[m] for m in members]
        unsorted_clusters.append(members)
    clustersizes = np.array([np.size(k,0) for k in unsorted_clusters])
    # sort by cluster sizes (largest cluster to smallest cluster)
    sorted_clusters = sorted(unsorted_clusters,key=len,reverse=True)
    
    # if there are the same pdbs in the clusters (from flipping), check for duplicates
    # also, only include if the cluster isn't a singleton
    unique_clusters = []
    for clus in sorted_clusters:
        #clus = list(k for k,_ in itertools.groupby(clus))
        labels = [str(k[0]) for k in clus] # label is the resn, ifgcount, vdmcount
        before = (len(labels)) # all labels
        after = (len(set(labels))) # unique label
        if after > 2:
            unique_clusters.append(clus)
            if before != after:
                print('There are duplicate pdbs in here from flipping coordinates')
                print('Example of something in this cluster:')
                print(labels[0])
    sorted_clusters = unique_clusters # rename

    if print_clusters == True:
        for i, members in enumerate(sorted_clusters):
            clus_num = str(i+1)
            for r in members:
                orig = r[1]
                name = '_'+r[1].split('/')[-1] 
                if BBorSC=='SC':
                    outfile = outputdir+'%s_cluster%s'%(vdmnames,clus_num+name)
                elif BBorSC=='N_CA' or BBorSC=='C_O':
                    outfile = outputdir+'cluster%s'%(clus_num+name)
                shutil.copy(orig, outfile)
    return sorted_clusters
    '''
