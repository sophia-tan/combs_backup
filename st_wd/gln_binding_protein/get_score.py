# commandline args: pdb of design, txtfile of vdms in this format:
#   resnum whitespace threeletterAA ifg 
#   Note: separate multiple ifgs with comma and NO WHITESPACE

import sys
from ScoringFunctions import *
from PyrosettaScores import *
import prody as pr
import freesasa, math
import numpy as np, pickle as pkl, pandas as pd

script, design_pdb, fg_vdm_txt= sys.argv

# 1) make pandas df to keep track of vdms and their scores
score_df = make_df(fg_vdm_txt)

# 2) add burial info
sasadict = cb_sasas(design_pdb, score_df)
score_df = scoring_sasa(sasadict, score_df)

# 3) add SS and rotamer info
score_df = pyrosetta_scores(design_pdb,score_df.index, score_df)
print(score_df)
