from kapp import RCAT
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import sys, os
sys.path.append(os.path.expanduser('~/git/global_functions'))
#from color import ColorMap
#from plot_types import cdf
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def boot_strap(x):
    x.dropna(inplace=True)
    median_values = np.zeros(1000)
    for i in xrange(1000):
        x = x.dropna()
        new_x = np.random.choice(x, len(x), replace=True)
        median_values[i] = np.median(new_x)
    return np.std(median_values)

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)

rc = RCAT(model)

rcat = rc.calculate_enzyme_rates()
rmax = rc.get_rcat_max(7)
kcat = rc.get_kcat_of_model_reactions()

reactions = rmax.index & kcat.index

rcat = rcat.loc[reactions]
rmax = rmax.loc[reactions]
kcat = kcat.loc[reactions]

gr = rc.growth_conditions.growth_rate_1_h
efficiency = rcat.div(kcat, axis=0).dropna()[gr.index].T
med = rcat.dropna(how='all').median()[gr.index]

colors = ['#324EAC' if rc.growth_conditions.growth_mode[c]=='chemostat' 
                            else '#3FBE57' for c in gr.index]

'''efficiency cdf plot'''
fig = plt.figure(figsize=(5,3.2))
ax = plt.axes()
#ax.set_axis_bgcolor((0.95,0.92,0.90))
medchemo = med[rc.growth_conditions[rc.growth_conditions.growth_mode=='chemostat'].index]
grchemo = gr[rc.growth_conditions[rc.growth_conditions.growth_mode=='chemostat'].index]
ax.scatter(grchemo, medchemo, zorder=10, c='#324EAC', edgecolor='none', s=50, alpha=0.75, label='chemostat')

medbatch = med[rc.growth_conditions[rc.growth_conditions.growth_mode=='batch'].index]
grbatch = gr[rc.growth_conditions[rc.growth_conditions.growth_mode=='batch'].index]
ax.scatter(grbatch, medbatch, zorder=10, c='#3FBE57', edgecolor='none', s=50, alpha=0.75, label='batch')

ax.axhline(kcat.median(), color='k', ls=':', lw=3, zorder=10)



for i, c in enumerate(gr.index):
    ax.errorbar(gr[c], med[c], boot_strap(rcat[c]), c=colors[i], zorder=10, alpha=0.5)

ax.text(.1, kcat.median()+0.25, 
         '$k_{\mathrm{cat}}$ median', 
         va='bottom', size=13)
ax.set_xticks(np.arange(0,0.8,0.2))
ax.set_yticks(np.arange(0,12.1,2))
ax.set_xticklabels(np.arange(0,0.8,0.2), size=12.5)
ax.set_yticklabels(np.arange(0,13,2), size=12.5)
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.set_xlabel(r'growth rate $[h^{-1}]$', size=13)
ax.set_ylabel('median $r_{\mathrm{cat}}$ $[s^{-1}]$', size=13)
plt.legend(loc=4, scatterpoints=1)
plt.tight_layout()
plt.savefig('../res/rcat_increase_with_gr.pdf')#