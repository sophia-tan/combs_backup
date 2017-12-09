import sys                              
from IFG_dicts import *
sys.path.append('/home/gpu/Sophia/combs/src/')
import combs
from sys import argv

script, ifg_dir, ifg = argv

ifg_sele = all_ifgs[ifg]
cb = combs.Comb(ifg_sele)

combs.make_bb_sc_rel_vdms(ifg_dir, cb)
