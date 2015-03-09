''' saturation and thermodynamics '''

from kapp import MM_KINETICS
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv 
import ast
from scipy.stats.mstats import gmean
import re

#import matplotlib.gridspec as gridspec
#plt.rcParams['text.usetex']=True
#plt.rcParams['text.latex.unicode']=True

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)
mm = MM_KINETICS(model)

KMs = pd.DataFrame.from_csv('../cache/KM_values.csv')
kcat = mm.get_kcat_of_model_reactions()
rcat = mm.calculate_enzyme_rates()
rmax = mm.get_rcat_max(7)
cmax = mm.get_cond_max(7)

relevant_reac = set(rmax.index & kcat.index & KMs.index)

i = 0

concentrations = mm.metab_conc
all_metabolites = mm.model_metabolites()

def calculate_saturation_effect(kms, stoicho, condition):
    
    metabolites = kms.index
    smul = 1
    for m in metabolites:
        cid = int(all_metabolites['kegg_id'][m][1:])
#        conc = 100 # uM
        if cid in concentrations[condition].dropna():
            conc = concentrations[condition][cid]*1000 
            s = (conc / kms[m])**stoicho[m]
            smul *= s
            return (smul) / (1 + smul)
    


carbon_cond = {r:mm.growth_conditions['carbon'][cmax[r]] for r in cmax.index}

dG, thermo = mm.backwards_reaction_effect()

S = pd.Series()
T = pd.Series()
i = 0
for r in model.reactions:
    if r.id in relevant_reac:
        stoicho = {k.name:-v for k,v in r.metabolites.iteritems() 
                                        if v<0 and k.name not in ['H2O', 'H+']}
        kms = KMs.loc[r.id].dropna()
        kms = kms[kms<0]*-1
        condition = carbon_cond[r.id]
        if condition in ['glc', 'ac', 'glyc']:
            if len(kms) == len(stoicho):
                S[r.id] = calculate_saturation_effect(kms, stoicho, condition)
            if r.id in thermo.index:
                T[r.id] = thermo.loc[r.id, condition]

T = T[T>-0.1]
T = np.abs(T)
S.dropna(inplace=True)
S = S.astype('float')

temp = T.index & S.index
S = S[temp]
T = T[temp]
ST = (S*T).dropna()

fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(8,3.5))
ax1.set_axis_bgcolor((0.95,0.92,0.90))
ax2.set_axis_bgcolor((0.95,0.92,0.90))
ax3.set_axis_bgcolor((0.95,0.92,0.90))

ax1.scatter(S, (rmax/kcat)[temp], s=50, alpha=0.4, edgecolor='none')            
ax1.set_xlim(1/1000.0, 3)
ax1.set_ylim(1/1000.0, 30)
ax1.axhline(1, 1e-3,3, ls=':', lw=3, color='r')
ax1.axvline(1, 0, 1e3, ls=':', lw=3, color='r')
ax1.set_ylabel(r'$r_{\rm{cat}}^{\rm{max}}$ / $k_{\rm{cat}}$', size=20)
ax1.set_xlabel('$ S $', size=20)
ax1.tick_params(axis='both', which='both', top='off', right='off')
ax1.plot(np.logspace(-3,0), np.logspace(-3,0), 'k:')

ax2.scatter(T, (rmax/kcat)[temp], s=50, color='y', alpha=0.4, edgecolor='none')            
ax2.set_xlim(1/1000.0, 3)
ax2.set_ylim(1/1000.0, 30)
ax2.axhline(1, 1e-3,3, ls=':', lw=3, color='r')
ax2.axvline(1, 0, 1e3, ls=':', lw=3, color='r')
ax2.set_xlabel('$ T $', size=20)
ax2.tick_params(axis='both', which='both', top='off', right='off')
ax2.plot(np.logspace(-3,0), np.logspace(-3,0), 'k:')

ax3.scatter(ST, (rmax/kcat)[ST.index], s=50, color='#00A352', alpha=0.4, edgecolor='none')            
ax3.set_xlim(1/1000.0, 3)
ax3.set_ylim(1/1000.0, 30)
ax3.axhline(1, 1e-3,3, ls=':', lw=3, color='r')
ax3.axvline(1, 0, 1e3, ls=':', lw=3, color='r')
ax3.set_xlabel(r'$ S $  x  $ T $', size=20)
ax3.tick_params(axis='both', which='both', top='off', right='off')
ax3.plot(np.logspace(-3,0), np.logspace(-3,0), 'k:')

r_to_b = mm.map_model_reaction_to_genes()

names = []
for r in ST.index:
    names.append(mm.gene_names[r_to_b[r]])

labels = list(ST[ST<0.4].index)  + ['PGI', 'PFK', 'AMAOTr', 'PTPATi']

for l in labels:
    if ST[l]:
        name = mm.gene_names[r_to_b[l]]
        if name == 'gltX':
            print name
        ax3.text(ST[l], (rmax/kcat)[l] , name,
                 ha='center', va='center', size=12)
        if S[l] < 0.4:
            if name != 'glyA':
                ax1.text(S[l], (rmax/kcat)[l] , name,
                         ha='center', va='center', size=12)
        if T[l] < 0.8:
            ax2.text(T[l], (rmax/kcat)[l] , name,
                     ha='center', va='center', size=12)
                 
ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xscale('log')
ax2.set_yscale('log')

ax3.set_xscale('log')
ax3.set_yscale('log')

plt.tight_layout()

plt.savefig('../res/saturation_and_thermodynamics.pdf')