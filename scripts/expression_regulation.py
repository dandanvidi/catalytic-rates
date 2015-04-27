from kapp import *
import sys, os
#sys.path.append(os.path.expanduser('~/git/component-contribution'))
sys.path.append(os.path.expanduser('~/git/cobrapy'))
from cobra.io.sbml import create_cobra_model_from_sbml_file
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from itertools import combinations
from scipy.optimize import curve_fit

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
rc = RCAT(model)

reactions = rc.calculate_enzyme_rates().index
gc = rc.growth_conditions
vilu = gc[gc.source=='vilu']
hein = gc[gc.source=='heinmann']
gm = gc.growth_mode
V = rc.V_data.loc[reactions]
V.replace([0.0, -0.0], np.nan, inplace=True)

E = rc.E_data.loc[reactions]

def draw_fold_change_hist(ax, array, cmap):
    ax.hist(array, histtype='stepfilled', 
            bins=np.arange(-1, np.max(array) + 0.1, 0.25),                                                         
            color='0.5', edgecolor='none', alpha=0.6)
    ax.axvline(array.median(), 0, 1, c=cmap(array.median()), zorder=10, ls='-', lw=4)
#    ax.axvspan(array.median(), 1, alpha=0.3, color='#cdb79e', zorder=3)
#    ax.axvline(1, 0, 1, c=cmap(0.99), zorder=10, ls='--', lw=3)
#    ax.axvline(0, 0, 1, c='k', zorder=10, ls='-', lw=1)
    
def get_flux_and_expression_fold_change(growthconditions, expression, fluxes):
    growthrates = growthconditions.growth_rate_1_h
    reactions = expression.index & fluxes.index
    combs = [(i,j) for (i,j) in combinations(growthrates.index, 2) 
                                    if growthrates[j] > growthrates[i]]
    mu_ratio = np.array([growthrates[j]/growthrates[i] for (i,j) in combs])
    mu2_to_mu1 = pd.DataFrame(index=combs, columns=['mu_ratio'])
    mu2_to_mu1.mu_ratio = mu_ratio

    E2_to_E1 = pd.DataFrame(index=reactions, columns=combs)
    V2_to_V1 = pd.DataFrame(index=reactions, columns=combs)    
    for i,j in combs:

        V2_to_V1[(i,j)] = np.log2(fluxes[j]/fluxes[i])
        E2_to_E1[(i,j)] = np.log2(expression[j]/expression[i])
        
    return E2_to_E1, V2_to_V1, mu2_to_mu1.T


fig = plt.figure(figsize=(6,6))
cm = plt.cm.get_cmap('Reds')
ax= plt.axes()
dE_vilu, dV_vilu, dmu_vilu = get_flux_and_expression_fold_change(vilu, E, V)
dE_hein, dV_hein, dmu_hein = get_flux_and_expression_fold_change(hein, E, V)

dE = dE_vilu.join(dE_hein)
dV = dV_vilu.join(dV_hein)
dmu = dmu_vilu.join(dmu_hein)

dE_to_dV = dE/dV
dE_to_dV.replace([np.inf, -np.inf], np.nan, inplace=True)
dE_to_dV.dropna(how='any', inplace=True)

popt, pcov = curve_fit(lambda x, a:a*x, dV.median().values, dE.median().values)

sc = plt.scatter(dV.median(), dE.median(), 
                 s=40, edgecolor='0.5',
                 cmap=cm ,c=dE_to_dV.median(), vmin=0, vmax=1) 

yerr = dE.std()/np.sqrt(dE.shape[0])
xerr = dV.std()/np.sqrt(dV.shape[0])

plt.errorbar(dV.median(), dE.median(), yerr=yerr, xerr=xerr, zorder=0, lw=1.2, c='0.5', fmt='o') 

axrng = [-0.5,2.5]
plt.plot(axrng,axrng, 'k:', lw=3)
#plt.plot(axrng,popt*axrng, c=cmap(dE_to_dV.median().median()), ls='-', lw=2, zorder=0)

#cb = plt.colorbar(orientation='horizontal')
plt.xlabel('$\Delta \log_2 (v)$ median', size=20)
plt.ylabel('$\Delta \log_2 (E)$ median', size=20)
ax.set_xticklabels(np.arange(-0.5,3,0.5), size=14)
ax.set_yticklabels(np.arange(-0.5,3,0.5), size=14)
plt.xlim(axrng)
plt.ylim(axrng)
ax.tick_params(axis='both', which='both', top='off', right='off')
plt.tight_layout()
plt.savefig('../res/expression_regulation_scatter.svg')


fig, axes = plt.subplots(2,1, figsize=(4.5,6), sharex=True, sharey=True)
for i, c in enumerate([('acetate_heinmann', 'glucose_heinmann'),
                       ('Chemostat_vilu_011', 'Chemostat_vilu_049')]):
    data = dE_to_dV[c]
    print data.mean()
    draw_fold_change_hist(axes[i], data, cmap)
        
axes[1].set_xticklabels(np.arange(-1,3,0.5), size=14)
axes[1].set_xlim(-1,2.5)

axes[0].set_yticks(np.arange(0,50,10))
axes[0].set_yticklabels(np.arange(0,50,10), size=14)
axes[1].set_yticklabels(np.arange(0,50,10), size=14)

axes[0].tick_params(axis='both', which='both', top='off', right='off')
axes[1].tick_params(axis='both', which='both', top='off', right='off')

axes[0].set_ylabel('reactions', size = 20)
axes[1].set_ylabel('reactions', size = 20)

axes[1].set_xlabel('$\Delta \log (E) / \Delta \log (v)$', size = 20)
plt.tight_layout()
plt.savefig('../res/expression_regulation_hists.svg')