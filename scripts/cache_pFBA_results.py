from kapp import RCAT
import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
'''
    GENERATE FLUXOMICS CACHE:
    
    preforms FBA with minimization of the total sum of fluxes
    and export the results to a csv file.
    
    Arguments:
        growth conditions
    
    Retunrs:
        csv file with flux distributions across conditions 
        in molecules_per_second_per_gCD

    While pFBA is calculated across the entire metabolic network, 
    the output csv file does not include reactions which can be catalyzed
    by nore then a single gene, to allow a one-to-one relationship between
    genes, fluxes and enzyme abundances    
'''

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)

reactions = [r.id for r in model.reactions]
gc = RCAT(model).growth_conditions.T.to_dict()

fluxes = pd.DataFrame(index=reactions, columns=RCAT(model).growth_conditions.index)

for i, (c, v) in enumerate(gc.iteritems()):
    
    model = create_cobra_model_from_sbml_file(model_fname)
    rate = RCAT(model)
    fluxes[c] = rate._calculate_pFBA(title=c, growth_params=v)

#remove all zero rows
#fluxes = fluxes.loc[(fluxes!=0).any(axis=1)]

## add reaction names and gene names
reactions = pd.DataFrame(RCAT(model).map_model_reaction_to_genes())
reactions.rename(columns={0: 'bnumber'}, inplace=True) 
reactions.index.name = 'reaction'

names = pd.DataFrame.from_dict(RCAT(model).gene_names, orient='index')
names.name='gene_name'

reactions = reactions.join(names, on='bnumber') 

fluxes = reactions.join(fluxes, how='left')
fluxes.index.name = 'reaction'
fluxes.rename(columns={0: 'gene_name'}, inplace=True) 

# export to csv
fluxes.to_csv('cache/fluxes_[molecules_per_second_per_gCDW].csv')
