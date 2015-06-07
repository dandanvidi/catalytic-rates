import pandas as pd
import numpy as np
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from copy import deepcopy
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

growth_conditions_fname = "../data/new_dataset/growth_conditions.csv"
gc = pd.DataFrame.from_csv(growth_conditions_fname)

##############################################
'''prepare proteomics for growth conditions'''        
##############################################

heinmann_fname = "../data/new_dataset/heinmann_copies_per_cell.csv"
heinmann = pd.read_csv(heinmann_fname) # copies/cell
uni_to_b = {row[48:54]:row[0:5].split(';')[0].strip() 
                            for row in open("../data/all_ecoli_genes.txt", 'r')}
heinmann.replace(to_replace={'uniprot':uni_to_b}, inplace=True)
heinmann.set_index('uniprot', inplace=True)                                
heinmann.index.name = 'bnumber'

cell_volume = gc['Single cell volume [fl]'] # cell volume in fl
cell_volume = cell_volume / 2 # the data in the table of cell volumes is biased by a factor of 2

ch = list(gc[gc.source=='heinmann'].index)

heinmann[ch] = heinmann[ch] / cell_volume[ch] # copies/fl

vilu_fname = "../data/new_dataset/valgepea_copies_per_fl.csv"
vilu = pd.DataFrame.from_csv(vilu_fname) # copies/fl[cell]

peebo_fname = "../data/new_dataset/peebo_copies_per_fl.csv"
peebo = pd.DataFrame.from_csv(peebo_fname) # copies/fl[cell]

''' now all data sets are in copies/fl cytoplasm so we can merge them'''
proteomics = heinmann.join(vilu, how='outer').astype('float')
proteomics = proteomics.join(peebo, how='outer').astype('float')
proteomics.dropna(how='all', inplace=True) # remove "all NaN" rows

''' export results '''
proteomics = proteomics[gc.index]
proteomics.to_csv('../data/proteomics[copies_fl].csv')

########################################
'''prepare pFBA for growth conditions'''        
########################################
        
def perform_pFBA(model, cs, gr, ur):

    model = deepcopy(model)
    convert_to_irreversible(model)            

    rxns = dict([(r.id, r) for r in model.reactions])
    rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
    try:
        rxns['EX_' + cs + '_e'].lower_bound = -ur # redefine sole carbon source uptake reaction in mmol/gr/h
    except:
        print cs, ur
        rxns['EX_glc_e'].lower_bound = -ur
    rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = gr            
    print "solving pFBA",
    optimize_minimal_flux(model, already_irreversible=True)
    
    flux_dist = pd.DataFrame(model.solution.x_dict.items()).set_index(0)
    
    return flux_dist    

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)
reactions = map(lambda x: x.id, model.reactions)
fluxes = pd.DataFrame(index=reactions, columns=gc.index)

for c in gc.iterrows():
    model = create_cobra_model_from_sbml_file(model_fname)
    convert_to_irreversible(model)
    gr = c[1]['growth rate (h-1)']
    cs = c[1]['carbon source']
    ur = c[1]['uptake rate [mmol/gCDW/h]']
    if np.isnan(ur):
        ur = 18.5
    model = create_cobra_model_from_sbml_file(model_fname)
    fluxes[c[0]] = perform_pFBA(model, cs, gr, ur)
    print "- %s" %c[0]

fluxes.index.name = 'reaction'
''' export results '''
fluxes.to_csv('../data/flux[mmol_h_gCDW].csv')
