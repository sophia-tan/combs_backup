import sys
sys.path.append('/Users/npolizzi/Projects/combs/src/')
import combs
import prody as pr
import pickle
# import collections
import networkx as nx
# import time


def pickle_load(file_path):
    with open(file_path, "rb") as f:
        return pickle.load(f)

indir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170915/var1/'
# relvdm_carboxamide = pickle_load(indir + 'relvdm_carboxamide.pickle')
# relvdm_carboxylate = pickle_load(indir + 'relvdm_carboxylate.pickle')
# relvdm_amino = pickle_load(indir + 'relvdm_amino.pickle')
print('loading sample pickle')
sample = pickle_load(indir + 'sample.pickle')
print('making pose metrics')
sample.pose_metrics()
# pose = sample.poses[31]
i = 1
for ind in sample.ranked_poses_by_score[:4]:
    sample.poses[ind].supergraph = nx.Graph()
    sample.poses[ind].site_graph = {}
    sample.poses[ind].sub_poses = []
    print('making subposes of pose ' + str(ind))
    sample.poses[ind].make_sub_poses()
    outdir = '/Users/npolizzi/Projects/combs/results/glutamine_binding_protein/20170917/poses_conf0/pose' + str(i) + '/'
    sample.poses[ind].print_sub_poses(sample, outdir)
    pr.writePDB(outdir + 'lig.pdb', sample.poses[ind].ligand)
    i += 1










# relvdm_carboxamide.hotspot_subgraphs = list(
#     sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxamide.hotspot_graph) if len(comp) > 4),
#            key=len, reverse=True))
#
# relvdm_carboxylate.hotspot_subgraphs = list(
#     sorted((comp for comp in nx.connected_component_subgraphs(relvdm_carboxylate.hotspot_graph) if len(comp) > 4),
#            key=len, reverse=True))
#
# relvdm_amino.hotspot_subgraphs = list(
#     sorted((comp for comp in nx.connected_component_subgraphs(relvdm_amino.hotspot_graph) if len(comp) > 4),
#            key=len, reverse=True))

# sample.ligand_conformers = [pr.parsePDB('/Users/npolizzi/Desktop/gln227_mm_good.pdb')]
# sample.set_ligand_ifg_correspondence(relvdm_carboxamide, 'CD CD NE2', 'CD NE2 NE2')
# sample.set_ligand_ifg_correspondence(relvdm_carboxylate, 'C C OXT', 'CD OE2 OE2')
# sample.set_ligand_ifg_correspondence(relvdm_amino, 'CA CA N', 'CE NZ NZ')
# sample.set_ligand_ifg_coords(relvdm_carboxamide)
# sample.set_ligand_ifg_coords(relvdm_carboxylate)
# sample.set_ligand_ifg_coords(relvdm_amino)
# sample.set_rel_vdms([relvdm_carboxamide, relvdm_amino, relvdm_carboxylate])
# sample.find_ligand_cst_hotspots(tol=0.2)

# sample.set_protein_bb_coords(clash_cutoff=2.5)
# sample.make_densities()
# sample.set_lig_ifg_dict(relvdm_carboxamide, 'CG CD OE1 NE2')
# sample.set_lig_ifg_dict(relvdm_carboxylate, 'CA C O OXT')
# sample.set_lig_ifg_dict(relvdm_amino, 'C CA N')  # NOTE: could have chosen C CA N here
# sample.set_lig_coords()
# sample.fit_lig_to_density()





