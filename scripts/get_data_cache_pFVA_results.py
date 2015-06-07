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


from kapp import RCAT
import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)

reactions = [r.id for r in model.reactions]


gc = RCAT(model).growth_conditions.T.to_dict()

fmin = pd.DataFrame(index=reactions, 
                    columns=list(RCAT(model).growth_conditions.index))
fmax = pd.DataFrame(index=reactions, columns=RCAT(model).growth_conditions.index).astype('float')

fmin.index.name = 'reactions'
fmax.index.name = 'reactions'

for i, (c, v) in enumerate(gc.iteritems()):
    
    model = create_cobra_model_from_sbml_file(model_fname)
    rate = RCAT(model)
    fva = rate.calculate_pFVA(title=c, growth_params=v)
    fmin[c] = [fva[r]['minimum'] for r in reactions]
    fmax[c] = [fva[r]['maximum'] for r in reactions]

# convert units from mmol/gCDW/h to molecules/gCDW/sec
fmin = fmin * 6.02214129e23 / 1000 / 3600 
fmax = fmax * 6.02214129e23 / 1000 / 3600 

# add bnumbers and gene names to dataframe  
bnumbers = RCAT(model).map_model_reaction_to_genes()
bnumbers.name='bnumber'

gene_names = pd.Series(RCAT(model).gene_names)
gene_names.name='gene_name'

genes = pd.DataFrame(bnumbers).join(gene_names, on='bnumber')

fmin = fmin.join(genes, how='left')
fmax = fmax.join(genes, how='left')

fmin.set_index('bnumber', append=True, inplace=True)
fmax.set_index('bnumber', append=True, inplace=True)

fmin.set_index('gene_name', append=True, inplace=True)
fmax.set_index('gene_name', append=True, inplace=True)

fmin.to_csv('cache/fmin_[molecules_per_second_per_gCDW].csv')
fmax.to_csv('cache/fmax_[molecules_per_second_per_gCDW].csv')
