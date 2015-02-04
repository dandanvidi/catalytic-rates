from kapp import RCAT
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def boot_strap(x):
    median_values = np.zeros(1000)
    for i in xrange(1000):
        x = x.dropna()
        new_x = np.random.choice(x, len(x), replace=True)
        median_values[i] = np.median(new_x)
    return np.std(median_values)

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)

rc = RCAT(model)

rcat = rc.calculate_enzyme_rates()
rmax = rc.get_rcat_max(7)
kcat = rc.get_kcat_of_model_reactions()
E = rc.E_data
gr = rc.growth_conditions.growth_rate_1_h
efficacy = rcat.div(kcat, axis=0).dropna().median()[gr.index]
med = rcat.dropna(how='all').median()[gr.index]

'''efficacy scatter plot'''
fig = plt.figure(figsize=(7,4.5))
ax = fig.add_subplot(111, axisbg='0.95')

ax.scatter(gr, efficacy, s=100, edgecolor='w', color='b', alpha=0.5)

for c in efficacy.index:
    plt.errorbar(gr[c], efficacy[c], color='b', alpha=0.5, 
                 yerr=boot_strap(rcat.div(kcat, axis=0).dropna()[c]))

max_efficacy = (rmax/kcat).dropna().median()
ax.hlines(max_efficacy, 0, 0.7, lw=3, color='r', linestyles=':')
ax.text(.1, max_efficacy+0.01, 
         r'$r_{\mathrm{cat}}^{\mathrm{max}}$  /  $k_{\mathrm{cat}}$', 
         va='bottom', size=15)

ax.set_xlabel(r'growth rate $\left[s^{-1}\right]$', size=15)
ax.set_ylabel('median catalytic efficacy \n [$r_{\mathrm{cat}} / k_{\mathrm{cat}}]$', size=15)

ax.tick_params(axis='x', which='both', top='off', bottom='on')
ax.tick_params(axis='y', which='both', left='on', right='off')

ax.set_xlim(0,0.7)
ax.set_ylim(0,1)

plt.tight_layout()
plt.savefig('res/catalytic_efficacy.pdf')

'''rcat increase scatter plot'''
fig = plt.figure(figsize=(7,4.5))
ax = fig.add_subplot(111, axisbg='0.95')

ax.scatter(gr, med, s=100, edgecolor='w', color='#9966FF')

for c in efficacy.index:
    plt.errorbar(gr[c], med[c], color='#9966FF', 
                 yerr=boot_strap(rcat.dropna(how='all')[c]))

ax.hlines(rmax.median(), 0, 0.7, lw=3, color='k', linestyles=':')
ax.text(.1, rmax.median()+0.01, 
         r'$r_{\mathrm{cat}}^{\mathrm{max}}$', 
         va='bottom', size=15)

ax.hlines(kcat.median(), 0, 0.7, lw=3, color='r', linestyles=':')
ax.text(.1, kcat.median()+0.01, 
         r'$k_{\mathrm{cat}}$', 
         va='bottom', size=15)

ax.set_xlabel(r'growth rate $\left[s^{-1}\right]$', size=15)
ax.set_ylabel('median catalytic rate $[s^{-1}]$', size=15)

ax.tick_params(axis='x', which='both', top='off', bottom='on')
ax.tick_params(axis='y', which='both', left='on', right='off')

ax.set_xlim(0,0.7)
#ax.set_ylim(0,1)

plt.tight_layout()
plt.savefig('res/rcat_increases_with_growth.pdf')