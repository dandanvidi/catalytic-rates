''' merge proteomics data sets '''

import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file

#model_fname = "data/iJO1366_curated.xml"
#model = create_cobra_model_from_sbml_file(model_fname)
#
h = pd.DataFrame.from_csv("data/Heinmann_proteomics.csv") # copies/cell
v = pd.DataFrame.from_csv("data/Vilu_proteomics.csv") # copies/fl[cell]

growth_conditions = pd.DataFrame.from_csv("data/growth_conditions.csv")

cell_volumes = growth_conditions.cell_volume_l # cell volume in liter
h = h*cell_volumes[h.columns]/1e-15 # convert from # copies/cell to copies/fl[cell]

p = h.join(v, how='outer')
p = p * 1e12 / 0.3 # convert from copies/fl(cytoplasme) to copies/gCDW
#
p.dropna(how='all', inplace=True)

