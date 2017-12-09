import sys                              
sys.path.append('/home/gpu/Sophia/combs/src/')
import combs

# stores input dicts for Combs objects
all_ifgs = {}
all_ifgs['carboxamide'] = {'ASN': 'CB CG OD1 ND2', 'GLN': 'CG CD OE1 NE2'}
all_ifgs['carboxylate'] = {'ASP': 'CB CG OD1 OD2', 'GLU': 'CG CD OE1 OE2'}
all_ifgs['amino'] = {'LYS': 'CE CD NZ'}

