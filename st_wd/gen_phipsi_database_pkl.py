from prody import *
import pickle as pkl
from sys import argv
import traceback
from SSFunctions import *
from PhiPsiFunctions import *
# need to provide textfile in command line!

# textfile is file that contains pdb and chains of all the pdbs in non-redundant PDB we're using
script, textfile = argv 

# compile lists of pdbs and chains
pdbnames = []
chains = []

with open(textfile) as fo:
    for line in fo:
        pdbnames.append(line[:4])
        chains.append(line[4])
zipped = zip(pdbnames, chains)

# get phi psi's for each pdb and then bin
counts_df = create_phi_psi_matrix()
for pdb, chain in zipped:
    try:
        print(pdb)
        phipsi = get_phi_psi_list(pdb, chain)
        for res in phipsi:
            counts_df = inc_phi_psi_df_for_db(counts_df, res[0], res[1])
    except Exception:
        print(traceback)
        print(pdb, 'pdb could not be parsed')

pkl.dump(counts_df, open('phipsi_bin_counts_database_df.pkl','wb'))
