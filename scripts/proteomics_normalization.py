from kapp import *
import sys, os
#sys.path.append(os.path.expanduser('~/git/component-contribution'))
sys.path.append(os.path.expanduser('~/git/cobrapy'))
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from itertools import combinations

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
rc = RCAT(model)

reactions = rc.calculate_enzyme_rates().index
gc = rc.growth_conditions
vilu = gc[gc.source=='vilu']
hein = gc[gc.source=='heinmann']
gm = gc.growth_mode
V = rc.V_data.loc[reactions]
V.replace([0.0, -0.0], np.nan, inplace=True)
E = rc.E_data
barlist = plt.bar(range(len(E.sum())), E.sum())

for i, c in enumerate(gc.index):
    if c in vilu.index:
        barlist[i].set_color('r')
