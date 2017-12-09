# vdmresidues must be comma separated with NO spaces.
# ex) 'ASP,GLU'
import sys                              
from IFG_dicts import *
sys.path.append('/home/gpu/Sophia/combs/src/')
import combs
from sys import argv
import pickle as pkl

script, combsoutputdir, ifg, vdmresidues= argv

combs.do_hier_clus_for_bb_and_sc(combsoutputdir, ifg, vdmresidues, print_clus=True)

