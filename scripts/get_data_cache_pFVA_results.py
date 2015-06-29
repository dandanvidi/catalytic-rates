from rcat import RCAT
import pandas as pd
import numpy as np
from cobra.core import Metabolite, Reaction
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from copy import deepcopy
from cobra.flux_analysis.variability import flux_variability_analysis

def perform_pFVA(model, cs, gr, ur, reactions, fraction_of_optimum=1.0):

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
    rxns['Ec_biomass_iJO1366_core_53p95M'].lower_bound = gr            

    fake = Metabolite(id='fake')
    model.add_metabolites(fake)        
            
    for r in model.reactions:
        r.add_metabolites({fake:1})
        
    flux_counter = Reaction(name='flux_counter')
    flux_counter.add_metabolites(metabolites={fake:-1})                

    model.add_reaction(flux_counter) 
    model.change_objective(flux_counter)
    
    print "solving pFVA"
    fva = flux_variability_analysis(model, 
                                    reaction_list=reactions, 
                                    objective_sense='minimize',
                                    fraction_of_optimum=fraction_of_optimum)
                                    
    return fva
#    flux_ranges = pd.Series(reactions)
#
#    for r, v in fva.iteritems():
#        flux_ranges[r] = v['maximum'] - v['minimum']
#
##        flux_ranges.to_csv("cache/reactions_to_pfva_ranges")
#
#    return flux_ranges


if __name__ == "__main__":

    R = RCAT()
    model_fname = "../data/iJO1366.xml"
    model = create_cobra_model_from_sbml_file(model_fname)
    convert_to_irreversible(model)
    fluxes = pd.DataFrame(index=R.rmaxn.index, columns=['minimum', 'maximum'])
    
    for c in R.gc.iterrows():
        reactions = R.rmaxn[R.rmaxn.condition==c[0]].index
        if len(reactions)!=0:
            model = create_cobra_model_from_sbml_file(model_fname)
            convert_to_irreversible(model)
            gr = c[1]['growth rate (h-1)']
            cs = c[1]['carbon source']
            ur = c[1]['uptake rate [mmol/gCDW/h]']
            if np.isnan(ur):
                ur = 18.5
            model = create_cobra_model_from_sbml_file(model_fname)
            fva = perform_pFVA(model, cs, gr, ur, reactions)
            for k, v in fva.iteritems():
                fluxes['minimum'][k] = v['minimum']
                fluxes['maximum'][k] = v['maximum']
    
    fluxes.to_csv('../cache/flux_variability[mmol_h_gCDW]_relaxation=0.csv')

    model = create_cobra_model_from_sbml_file(model_fname)
    convert_to_irreversible(model)
    fluxes = pd.DataFrame(index=R.rmaxn.index, columns=['minimum', 'maximum'])
    
    for c in R.gc.iterrows():
        reactions = R.rmaxn[R.rmaxn.condition==c[0]].index
        if len(reactions)!=0:
            model = create_cobra_model_from_sbml_file(model_fname)
            convert_to_irreversible(model)
            gr = c[1]['growth rate (h-1)']
            cs = c[1]['carbon source']
            ur = c[1]['uptake rate [mmol/gCDW/h]']
            if np.isnan(ur):
                ur = 18.5
            model = create_cobra_model_from_sbml_file(model_fname)
            fva = perform_pFVA(model, cs, gr, ur, reactions)
            for k, v in fva.iteritems():
                fluxes['minimum'][k] = v['minimum']
                fluxes['maximum'][k] = v['maximum']
    
    fluxes.to_csv('../cache/flux_variability[mmol_h_gCDW]_relaxation=0_1.csv')
    