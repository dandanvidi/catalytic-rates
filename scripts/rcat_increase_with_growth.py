from kapp import RCAT
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import sys, os
sys.path.append(os.path.expanduser('~/git/global_functions'))
from color import ColorMap
from plot_types import cdf
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

colors = ['#996633' if rc.growth_conditions.carbon[c]!='glc' else '#CC0066' for c in gr.index]

'''efficiency cdf plot'''
fig = plt.figure(figsize=(5,3.2))
ax = plt.axes()
ax.set_axis_bgcolor((0.95,0.92,0.90))
ax.scatter(gr, med, zorder=10, c=colors, edgecolor='none', s=100, alpha=0.75)
ax.axhline(kcat.median(), color='k', ls='-', lw=2, zorder=10)
for i, c in enumerate(gr.index):
    ax.errorbar(gr[c], med[c], boot_strap(rcat[c]), c=colors[i], zorder=10)
ax.grid(color='w', ls='-', lw=1.2, zorder=0)
ax.tick_params(color='k')
ax.set_xticks(np.arange(0,0.8,0.2))
ax.set_yticks(np.arange(0,11,2))
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.set_xlabel(r'growth rate $\left[h^{-1}\right]$', size=13)
ax.set_ylabel('$r_{\mathrm{cat}} / k_{\mathrm{cat}}$', size=18)
plt.tight_layout()
plt.savefig('../res/rcat_increase_with_gr.pdf')#
#ax.tick_params(axis='x', which='both', top='off', bottom='on')
#ax.tick_params(axis='y', which='both', left='on', right='off')
#
#ax.set_xlim(0,0.7)
#ax.set_ylim(0,1)
#
#plt.tight_layout()
#plt.savefig('../res/catalytic_efficacy.pdf')
#
#'''rcat increase scatter plot'''
#fig = plt.figure(figsize=(7,4.5))
#ax = fig.add_subplot(111, axisbg='0.95')
#
#ax.scatter(gr, med, s=100, edgecolor='w', color='#9966FF')
#
#for c in efficacy.index:
#    plt.errorbar(gr[c], med[c], color='#9966FF', 
#                 yerr=boot_strap(rcat.dropna(how='all')[c]))
#
#ax.hlines(rmax.median(), 0, 0.7, lw=3, color='k', linestyles=':')
#ax.text(.1, rmax.median()+0.01, 
#         r'$r_{\mathrm{cat}}^{\mathrm{max}}\/median\/across\/all\/enzymes}$', 
#         va='bottom', size=15)
#
#ax.hlines(kcat.median(), 0, 0.7, lw=3, color='r', linestyles=':')
#ax.text(.1, kcat.median()+0.01, 
#         r'$k_{\mathrm{cat}}$', 
#         va='bottom', size=15)
#
#ax.set_xlabel(r'growth rate $\left[h^{-1}\right]$', size=15)
#ax.set_ylabel('median catalytic rate $[s^{-1}]$', size=15)
#
#ax.tick_params(axis='x', which='both', top='off', bottom='on')
#ax.tick_params(axis='y', which='both', left='on', right='off')
#
#ax.set_xlim(0,0.7)
##ax.set_ylim(0,1)
#
#plt.tight_layout()
#plt.savefig('../res/rcat_increases_with_growth.pdf')