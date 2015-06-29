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
    
    rmax = R.rmax['rmax [s^-1]'][R.rmaxn.index]
    rmax_mz = calculate_new_rmax(main_zero)[R.rmaxn.index]
    rmax_mh = calculate_new_rmax(main_high)[R.rmaxn.index]
    
    fig = plt.figure(figsize=(6,4))
    ax = plt.axes()
    plt.axhline(y=1, c='0.5', lw=2, zorder=0) 
    plt.legend(loc=4, scatterpoints=1)
    
    ax.set_xlim(0, len(rmax))
    ax.set_ylim(0, 2)
    ax.set_ylabel('$r_\mathrm{max}$ ratio', size=20)
    ax.set_xlabel('enzyme-reaction pair', size=20)
    

#    y = rmax[rmax_ol.index]/rmax_ol    
#    x = pd.Series(index=y.index, data = range(len(rmax_ol)))
    
    labels = [r for r in rmax.index if abs(np.log2(rmax/rmax_mz)[r])>0.05]
    labels += [r for r in rmax.index if abs(np.log2(rmax/rmax_mh)[r])>0.05]
    labels = list(set(labels))
    for r in labels:
        x = list(rmax.index).index(r)
        y = (rmax_mh/rmax)[r]
        y1 = (rmax_mz/rmax)[r]
        plt.scatter(x, y, c='g', edgecolor='none', s=15, 
                label='zero ATP maintenance')
        plt.scatter(x, y1, c='r', edgecolor='none', s=15, 
                label='200% ATP maintenance')
        if y>1.05 or y<0.95:
            ax.text(x, y+0.05, R.map_reactions_to_gene_names()[r], va='bottom', ha='center', fontsize=14)
        elif y1>1.05 or y1<0.95:
            ax.text(x, y1+0.05, R.map_reactions_to_gene_names()[r], va='bottom', ha='center', fontsize=14)
            
        plt.plot([x,x], [min(y, y1),max(y, y1)], c='k')
                
    ax.tick_params(axis='both', which='both', top='off', right='off')   
    
    [tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
    plt.tight_layout()
    plt.savefig('%s/svg/sensitivity_analysis_for_FBA.svg'%R.path)   