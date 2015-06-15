import pandas as pd
import numpy as np
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from copy import deepcopy
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from get_data_flux_and_proteomics import perform_pFBA
from rcat import RCAT
import matplotlib.pyplot as plt
from figure_correlation_plot import add_labels
R = RCAT()

def ATP_maintenance():
    R.rxns['ATPM'].lower_bound = 0
    fluxes = pd.DataFrame(index=R.rxns.keys(), columns=R.gc.index)
    for c in R.gc.iterrows():
        gr = c[1]['growth rate (h-1)']
        cs = c[1]['carbon source']
        ur = c[1]['uptake rate [mmol/gCDW/h]']
        if np.isnan(ur):
            ur = 18.5
        fluxes[c[0]] = perform_pFBA(R.model, cs, gr, ur)
        print "- %s" %c[0]
    
    fluxes.index.name = 'reaction'
    ''' export results '''
    fluxes.to_csv('../cache/flux[mmol_h_gCDW]_maintenance_zero.csv')
    
    R.rxns['ATPM'].lower_bound = 3.15*2
    for c in R.gc.iterrows():
        gr = c[1]['growth rate (h-1)']
        cs = c[1]['carbon source']
        ur = c[1]['uptake rate [mmol/gCDW/h]']
        if np.isnan(ur):
            ur = 18.5
        fluxes[c[0]] = perform_pFBA(R.model, cs, gr, ur)
        print "- %s" %c[0]
    
    fluxes.index.name = 'reaction'
    fluxes.to_csv('../cache/flux[mmol_h_gCDW]_maintenance_high.csv')
    
    R.rxns['ATPM'].lower_bound = 3.15

def oxygen_uptake():
    fluxes = pd.DataFrame(index=R.rxns.keys(), columns=R.gc.index)
    for c in R.gc.iterrows():
    
        gr = c[1]['growth rate (h-1)']
        cs = c[1]['carbon source']
        ur = c[1]['uptake rate [mmol/gCDW/h]']
    
        oxy = R.v[c[0]]['EX_o2_e_reverse']
        R.rxns['EX_o2_e_reverse'].upper_bound = oxy*0.5
        if np.isnan(ur):
            ur = 18.5
        fluxes[c[0]] = perform_pFBA(R.model, cs, gr, ur)
        print "- %s" %c[0]
    #    
    fluxes.index.name = 'reaction'
    ''' export results '''
    fluxes.to_csv('../cache/flux[mmol_h_gCDW]_oxygen_low.csv')
    
    R.rxns['EX_o2_e_reverse'].lower_bound = 0
    R.rxns['EX_o2_e_reverse'].upper_bound = 1000

def calculate_new_rmax(fluxes):
    temp = R.v
    R.v = fluxes
    rcat = R.calculate_catalytic_rates()
    rmax = rcat.dropna(thresh=5).max(axis=1)
    R.v = temp
    return rmax

if __name__ == "__main__":
    
    main_zero = pd.DataFrame.from_csv('../cache/flux[mmol_h_gCDW]_maintenance_zero.csv')       
    main_high = pd.DataFrame.from_csv('../cache/flux[mmol_h_gCDW]_maintenance_high.csv')   
    low_oxy = pd.DataFrame.from_csv('../cache/flux[mmol_h_gCDW]_oxygen_low.csv')    
    
    rmax = R.rmax['rmax [s^-1]']    
    
    rmax_mz = calculate_new_rmax(main_zero)
    rmax_mh = calculate_new_rmax(main_high)
    rmax_ol = calculate_new_rmax(low_oxy)
    
    fig = plt.figure(figsize=(6,4))
    ax = plt.axes()
    plt.scatter(range(len(rmax)), rmax/rmax_mz, c='r', edgecolor='none', s=10, alpha=0.5, label='zero ATP maintenance')
    plt.scatter(range(len(rmax)), rmax/rmax_mh[rmax.index], c='g', edgecolor='none', alpha=0.5, s=10, label='200% ATP maintenance')
    plt.scatter(range(len(rmax_ol)), rmax[rmax_ol.index]/rmax_ol, c='b', edgecolor='none', alpha=0.5, s=10, label='50% oxygen uptake')
    ax.set_yscale('log', basey=2)
#    ax.set_yscale('log')
    plt.legend(loc=4, scatterpoints=1)
    
    ax.set_xlim(0, 250)
    ax.set_ylim(2**-6, 2**3)
    ax.set_ylabel('$r_\mathrm{max}$ ratio', size=15)
    ax.set_xlabel('enzyme-reaction pair', size=15)
    

    y = rmax[rmax_ol.index]/rmax_ol    
    x = pd.Series(index=y.index, data = range(len(rmax_ol)))
    
    labels = {k:v for k,v in R.map_reactions_to_gene_names().iteritems() if k in y.index
                and abs(np.log2(y[k])) > 1}
#    labels.pop("GART", None)
    add_labels(x, y, labels, ax, fig, hide_overlap=False)
#    ax.text(y.index['GART'], y['GART'], 'purT')
    
    plt.tight_layout()
    plt.savefig('../res/sensitivity_analysis_for_FBA.svg')   