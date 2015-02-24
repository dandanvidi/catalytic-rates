''' map proteomics data sets '''

from kapp import RCAT
import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
#
hein = pd.DataFrame.from_csv("data/Heinmann_proteomics.csv") # copies/cell
vilu = pd.DataFrame.from_csv("data/Vilu_proteomics.csv") # copies/fl[cell]

growth_conditions = pd.DataFrame.from_csv("data/growth_conditions.csv")

cell_volumes = growth_conditions.cell_volume_l # cell volume in liter
correction = 2 # the data in the table of cell volumes is biased by a factor of 2
hein = hein / ((cell_volumes[hein.columns]/correction) * 1e15) # convert from # copies/cell to copies/fl[cell]

expression = hein.join(vilu, how='outer')
expression = expression[growth_conditions.index]

''' now both data sets are in copies/fl cytoplasm '''
expression = (expression[expression>10])# remove genes with less then 10 copies per fl cytoplasm

rho = 1100 # average cell density gr/liter
DW_fraction = 0.3 # fraction of DW of cells
expression = expression * 1e15 / (rho * DW_fraction) # convert from copies/fl(cytoplasme) to copies/gCDW

expression.dropna(how='all', inplace=True)

## add reaction names and gene names
reactions = pd.DataFrame(RCAT(model).map_model_reaction_to_genes())
reactions.rename(columns={0: 'bnumber'}, inplace=True) 
reactions.index.name = 'reaction'

names = pd.DataFrame.from_dict(RCAT(model).gene_names, orient='index')
reactions.rename(columns={0: 'gene_name'}, inplace=True) 
names.name='gene_name'

reactions = reactions.join(names, on='bnumber') 

expression = reactions.join(expression, on='bnumber', how='left')
expression.rename(columns={0: 'gene_name'}, inplace=True) 

#expression.dropna(how='all', inplace=True)

# sort by growth rate
expression = expression[['bnumber', 'gene_name'] + 
                            list(RCAT(model).growth_conditions.index)]

# export to csv
expression.to_csv('cache/expression_[copies_per_gCDW].csv')
