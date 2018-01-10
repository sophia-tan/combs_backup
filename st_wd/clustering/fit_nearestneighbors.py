from os import listdir 
from os.path import isfile, join
from sys import argv
import pickle as pkl
import math
from sklearn.neighbors import NearestNeighbors

script, ifg = argv

num_atoms = {}
num_atoms['backboneCO'] = 2
num_atoms['serCOH'] = 2
num_atoms['thrCOH'] = 2
num_atoms['tyrCOH'] = 2
num_atoms['indole'] = 3
num_atoms['imidazole'] = 5
num_atoms['amino'] = 3
num_atoms['carboxamide'] = 4
num_atoms['carboxylate'] = 4
for dist in [1, 1.2, 1.4, 1.6, 1.8]:
    for typ in ['SC', 'N_CA', 'C_O']:
        path = '/home/gpu/Sophia/STcombs/20171118/%s/clusters/%s/pickle/'%(ifg,typ)
        onlyfiles = [f for f in listdir(path) if isfile(join(path,f))]
        for pklfile in onlyfiles:
            pklf = pkl.load(open(path+pklfile, 'rb'))
            ifgcoords = pklf[:,6]
            flat = [x.flatten() for x in ifgcoords]
            #dist = 0.7*math.sqrt(num_atoms[ifg])
            ##print(dist)
            nbrs = NearestNeighbors(n_neighbors=1,metric='euclidean',radius=dist)
            vdmname = pklfile.split('_')[0]
            x = nbrs.fit(flat)
            pkl.dump(x, open('NNfit_%s_with_%s_%s_%s.pkl' %(ifg,vdmname,typ,dist), 'wb'))
