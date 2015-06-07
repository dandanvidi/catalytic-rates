from rcat import RCAT
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from itertools import combinations
from scipy.optimize import curve_fit

R = RCAT()

reactions = R.rcat.index

valg = R.gc[R.gc.source=='valgepea']
pebo = R.gc[R.gc.source=='peebo']
hein = R.gc[R.gc.source=='heinmann']

v = R.convert_mmol_gCDW_h_to_molecules_gCDW_s(R.v).loc[reactions]
v.replace(0.0, np.nan, inplace=True)

p = v/R.rcat

def draw_foldchange_hist(ax, array, cmap):
    ax.hist(array, histtype='stepfilled', 
            bins=np.arange(-1, np.max(array) + 0.1, 0.25),                                                         
            color='0.5', edgecolor='none', alpha=0.6)
    ax.axvline(array.median(), 0, 1, c=cmap(array.median()), zorder=10, ls='-', lw=4)
    
def get_v_and_p_foldchange(v, p, conds):
    conds = conds.index
    v = v[conds]
    p = p[conds]

    gr = R.gc['growth rate (h-1)']
    combs = [(i,j) for (i,j) in combinations(conds, 2) if gr[j] > gr[i]]
    combs_all = [(i,j) for (i,j) in combinations(R.gc.index, 2) if gr[j] > gr[i]]    
    delta_mu = pd.Series(data = map(lambda x: np.log2(gr[x[1]]/gr[x[0]]), combs_all),
                         index = combs_all)
    delta_p = pd.DataFrame(index=reactions, columns=combs)
    delta_v = pd.DataFrame(index=reactions, columns=combs)
    for (i, j) in combs:
        delta_p[(i,j)] = np.log2(p[j] / p[i])
        delta_v[(i,j)] = np.log2(v[j] / v[i])
    return delta_p, delta_v, delta_mu
    

fig = plt.figure(figsize=(6,6))
cm = plt.cm.get_cmap('Reds')
ax= plt.axes()
dp_valg, dv_valg, dmu_valg = get_v_and_p_foldchange(v, p, valg)
dp_hein, dv_hein, dmu_hein = get_v_and_p_foldchange(v, p, hein)
dp_pebo, dv_pebo, dmu_pebo = get_v_and_p_foldchange(v, p, pebo)

dp = dp_valg.join(dp_hein).join(dp_pebo)
dv = dv_valg.join(dv_hein).join(dv_pebo)
dmu = dmu_pebo.loc[dv.columns]

dp_to_dv = (dp/dv).dropna(how='any')
#popt, pcov = curve_fit(lambda x, a:a*x, dV.median().values, dE.median().values)
#
sc = plt.scatter(dv.median(), dp.median(), 
                 s=40, edgecolor='0.5',
                 cmap=cm ,c=dp_to_dv.median(), vmin=0, vmax=1) 
#
#yerr = dE.std()/np.sqrt(dE.shape[0])
#xerr = dV.std()/np.sqrt(dV.shape[0])
#
#plt.errorbar(dV.median(), dE.median(), yerr=yerr, xerr=xerr, zorder=0, lw=1.2, c='0.5', fmt='o') 
#
axrng = [-0.5,2.5]
plt.plot(axrng,axrng, 'k:', lw=3)
#plt.plot(axrng,popt*axrng, c=cmap(dE_to_dV.median().median()), ls='-', lw=2, zorder=0)
#
##cb = plt.colorbar(orientation='horizontal')
#plt.xlabel('$\Delta \log_2 (v)$ median', size=20)
#plt.ylabel('$\Delta \log_2 (E)$ median', size=20)
#ax.set_xticklabels(np.arange(-0.5,3,0.5), size=14)
#ax.set_yticklabels(np.arange(-0.5,3,0.5), size=14)
plt.xlim(axrng)
plt.ylim(axrng)
#ax.tick_params(axis='both', which='both', top='off', right='off')
#plt.tight_layout()
#plt.savefig('../res/expression_regulation_scatter.svg')


#fig, axes = plt.subplots(2,1, figsize=(4.5,6), sharex=True, sharey=True)
#for i, c in enumerate([('acetate_heinmann', 'glucose_heinmann'),
#                       ('Chemostat_vilu_011', 'Chemostat_vilu_049')]):
#    data = dE_to_dV[c]
#    print data.mean()
#    draw_fold_change_hist(axes[i], data, cmap)
#        
#axes[1].set_xticklabels(np.arange(-1,3,0.5), size=14)
#axes[1].set_xlim(-1,2.5)
#
#axes[0].set_yticks(np.arange(0,50,10))
#axes[0].set_yticklabels(np.arange(0,50,10), size=14)
#axes[1].set_yticklabels(np.arange(0,50,10), size=14)
#
#axes[0].tick_params(axis='both', which='both', top='off', right='off')
#axes[1].tick_params(axis='both', which='both', top='off', right='off')
#
#axes[0].set_ylabel('reactions', size = 20)
#axes[1].set_ylabel('reactions', size = 20)
#
#axes[1].set_xlabel('$\Delta \log (E) / \Delta \log (v)$', size = 20)
#plt.tight_layout()
#plt.savefig('../res/expression_regulation_hists.svg')