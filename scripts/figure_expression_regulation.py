from rcat import RCAT
import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
from itertools import combinations
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

def boot_strap(x):
    x.dropna(inplace=True)
    median_values = np.zeros(1000)
    for i in xrange(1000):
        x = x.dropna()
        new_x = np.random.choice(x, len(x), replace=True)
        median_values[i] = np.median(new_x)
    return np.std(median_values)
    
def draw_foldchange_hist(ax, array):
    ax.hist(array, histtype='stepfilled', 
            bins=np.arange(-1, np.max(array) + 0.1, 0.25),                                                         
            color='#CCCCB2', edgecolor='none')
    ax.axvline(array.median(), 0, 1, c='#438343', zorder=10, ls='-', lw=4)

def get_foldchange(v, p, conds):

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

R = RCAT()
fontsize = 20
reactions = R.rcat.index
gr = R.gc['growth rate (h-1)']

valg = R.gc[R.gc.source=='valgepea']
pebo = R.gc[R.gc.source=='peebo']
hein = R.gc[R.gc.source=='heinmann']

v = R.convert_mmol_gCDW_h_to_molecules_gCDW_s(R.v).loc[reactions]
v.replace(0.0, np.nan, inplace=True)

p = v/R.rcatn

dp_valg, dv_valg, dmu_valg = get_foldchange(v, p, valg.index)
dp_hein, dv_hein, dmu_hein = get_foldchange(v, p, hein.index)
dp_pebo, dv_pebo, dmu_pebo = get_foldchange(v, p, pebo.index)

dp = dp_valg.join(dp_hein).join(dp_pebo)
dv = dv_valg.join(dv_hein).join(dv_pebo)
dmu = dmu_pebo.loc[dv.columns]

dp_to_dv = (dp/dv).dropna(how='any')

plt.figure(figsize = (10,8.5))
cm = plt.cm.get_cmap('Reds')

gs = gridspec.GridSpec(3, 3)
ax1 = plt.subplot(gs[0:2, 0:2])
ax2 = plt.subplot(gs[0,2])
ax3 = plt.subplot(gs[1,2])

scat = R.convert_copies_gCDW_to_copies_fl(p)
ax1.scatter(gr, scat.median(), s=80, c='#65C565',
            edgecolor='none')

for c in gr.index:
    ax1.errorbar(gr[c], scat[c].median(), boot_strap(scat[c]), c='#65C565', zorder=0)
    

hist2 = ('Acetate', 'Glucose + salts')
hist3 = ('Chemostat u=0.11', 'Chemostat u=0.49')
draw_foldchange_hist(ax2, dp_to_dv[hist2])
draw_foldchange_hist(ax3, dp_to_dv[hist3])



ax1.set_ylim(0, 700)
#
[tick.label.set_fontsize(fontsize) for tick in ax1.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax1.yaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax2.yaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax3.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax3.yaxis.get_major_ticks()]
#
ax1.set_yticklabels(['0', '', '200', '', '400', '', '600', ''])
ax1.set_xticklabels(['', '', '0.2', '', '0.4', '', '0.6', ''])
#
ax1.set_xlabel(r'growth rate $[h^{-1}]$', size=fontsize)
ax1.set_ylabel(r'$E\,\left[\frac{\rm{copies}}{\rm{fL}}\right]$', size=fontsize)
ax1.tick_params(axis='both', which='both', top='off', right='off')

ax2.set_yticks(np.arange(5,45,10))
ax2.set_xticklabels([''])
ax2.set_xlim(-1.1, 2.1)
ax2.set_ylim(0, 25)
ax2.set_ylabel('reactions', size=fontsize)
ax2.tick_params(axis='both', which='both', top='off', right='off')
#ax2.set_yticklabels([''])
#
ax3.set_yticks(np.arange(5,45,10))

ax3.set_xlim(-1.1, 2.1)
ax3.set_xticks(np.arange(-1,2.1, 0.5))
ax3.set_xticklabels(['-1', '', '0', '', '1', '', '2'])


ax3.set_ylim(0, 25)
#
#

#

#
ax3.set_xlabel(r'$\Delta \log (E) / \Delta \log (v)$', size=fontsize)
ax3.set_ylabel('reactions', size=fontsize)



ax3.tick_params(axis='both', which='both', top='off', right='off')

plt.tight_layout()
plt.savefig('%s/svg/expression_regulation.svg'%R.path)














#
#
#fig = plt.figure(figsize=(6,6))
#ax = plt.axes()
#ax.scatter(gr, R.convert_copies_gCDW_to_copies_fl(p.median()))
#
#

#    
#def get_foldchange(v, p, conds):
#    v = v[conds]
#    p = p[conds]
#
#    gr = R.gc['growth rate (h-1)']
#    combs = [(i,j) for (i,j) in combinations(conds, 2) if gr[j] > gr[i]]
#    combs_all = [(i,j) for (i,j) in combinations(R.gc.index, 2) if gr[j] > gr[i]]    
#    delta_mu = pd.Series(data = map(lambda x: gr[x[1]]/gr[x[0]], combs_all),
#                         index = combs_all)
#    delta_p = pd.DataFrame(index=reactions, columns=combs)
#    delta_v = pd.DataFrame(index=reactions, columns=combs)
#    for (i, j) in combs:
#        delta_p[(i,j)] = p[j] / p[i]
#        delta_v[(i,j)] = v[j] / v[i]
#    return delta_p, delta_v, delta_mu
#    
#
#dp_valg, dv_valg, dmu_valg = get_foldchange(v, p, valg.index)
#dp_hein, dv_hein, dmu_hein = get_foldchange(v, p, hein.index)
#dp_pebo, dv_pebo, dmu_pebo = get_foldchange(v, p, pebo.index)

#fig = plt.figure(figsize=(7,7))
#
#ax= plt.axes()
#
#
#dp = dp_valg.join(dp_hein).join(dp_pebo)
#dv = dv_valg.join(dv_hein).join(dv_pebo)
#dmu = dmu_pebo.loc[dv.columns]
#
#dp_to_dv = (dp/dv).dropna(how='any')
##popt, pcov = curve_fit(lambda x, a:a*x, dV.median().values, dE.median().values)
##
#
#for (c1, c2) in dv.columns:
#    if R.gc['carbon source'][c1] != R.gc['carbon source'][c2]:
#
#        plt.scatter(dv.median()[(c1,c2)], dp.median()[(c1,c2)], 
#                 s=40, edgecolor='none', marker='o', alpha=0.5,
#                 vmin=0, vmax=1) 
#    else:
#        plt.scatter(dv.median()[(c1,c2)], dp.median()[(c1,c2)], 
#                    s=40, edgecolor='none', marker='o', alpha=0.5,
#                    vmin=0, vmax=1, color='r') 
#
##
#yerr = dp.std()/np.sqrt(dp.shape[0])
#xerr = dv.std()/np.sqrt(dv.shape[0])
##
##plt.errorbar(dv.median(), dp.median(), yerr=yerr, xerr=xerr, zorder=0, lw=1.2, c='0.5', fmt='o') 
##
#axrng = [-0.3,2.5]
#plt.plot(axrng,axrng, 'k:', lw=3)
##plt.plot(axrng,popt*axrng, c=cmap(dE_to_dV.median().median()), ls='-', lw=2, zorder=0)
##
###cb = plt.colorbar(orientation='horizontal')
#plt.xlabel('$\Delta \log_2 (v)$ median', size=fontsize)
#plt.ylabel('$\Delta \log_2 (E)$ median', size=fontsize)
##ax.set_xticklabels(np.arange(-0.5,3,0.5), size=14)
##ax.set_yticklabels(np.arange(-0.5,3,0.5), size=14)
#plt.xlim(axrng)
#plt.ylim(-.3, 1)
#ax.tick_params(axis='both', which='both', top='off', right='off')
#[tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
#[tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]
#
#plt.tight_layout()
#plt.savefig('%s/svg/expression_regulation_scatter.svg'%R.path)


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