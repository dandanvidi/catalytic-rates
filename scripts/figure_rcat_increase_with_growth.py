from rcat import RCAT
import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append(os.path.expanduser('~/git/across-projects'))
from plot_types import cdf
from color import ColorMap

fontsize = 20

def boot_strap(x):
    x.dropna(inplace=True)
    median_values = np.zeros(1000)
    for i in xrange(1000):
        x = x.dropna()
        new_x = np.random.choice(x, len(x), replace=True)
        median_values[i] = np.median(new_x)
    return np.std(median_values)


R = RCAT()

kcat = R.kcat['kcat [s^-1]']
x = R.gc['growth rate (h-1)']
rcatn = R.rcatn
y = rcatn.median()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5,6))

cm = plt.cm.get_cmap('Blues')
conds = ['Chemostat u=0.11', 'Chemostat u=0.20', 'Chemostat u=0.31', 'Chemostat u=0.51']
i = 0
for c in conds:
    m = y[c]
    label = '$\mu=%.01f\,h^{-1}$'%x[c]
    cdf(rcatn[c], color=cm(x[c]*1.5/x.max()), ax=ax1, label=label, lw=2.5)
    ax1.plot([m, m], [0, 0.5], c=cm(x[c]*1.5/x.max()), ls='-')
    props = dict(boxstyle='round', facecolor=cm(x[c]*0.8/x.max()))
    ax1.text(m, 0.025+i, '%0.1f'%m, bbox=props, ha='center', size=fontsize/1.7)    
    i += 0.07
    
cdf(kcat, color='k', ax=ax1, label=r'$k_{\mathrm{cat}}$', lw=2.5)
ax1.plot([kcat.median(), kcat.median()], [0, 0.5], c='k', ls='-')
props = dict(boxstyle='round', facecolor='0.7')
ax1.text(kcat.median(), 0.02+i, '%0.1f'%kcat.median(), bbox=props, ha='center', size=fontsize/1.7)    

ax1.set_xscale('log')
ax1.set_xlim(1e-2, 1e3)
ax1.set_xlabel(r'$r_{\mathrm{cat}}\,[s^{-1}]$', size=fontsize)
ax1.set_ylabel('cumulative distribution', size=fontsize)
ax1.legend(loc=0, fontsize=fontsize/1.1, frameon=False, handlelength=1)


ax2.set_xticklabels(['', '0.2', '', '0.4', '', '0.6', '', '0.8'])

for i, c in enumerate(x.index):
    if c not in conds:
        col = '1'
    else:
        col = cm(x[c]/x.max())
    ax2.scatter(x[c], y[c], zorder=10, c=col, 
            edgecolor='k', s=75)
    ax2.errorbar(x[c], y[c], boot_strap(rcatn[c]), c='k', zorder=0, alpha=0.5)
        
ax2.set_ylim(0,9)
#ax.set_ylim(0,0.5)
ax1.tick_params(axis='both', which='both', top='off', right='off')
ax2.tick_params(axis='both', which='both', top='off', right='off')

[tick.label.set_fontsize(fontsize) for tick in ax1.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax1.yaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax2.xaxis.get_major_ticks()]
[tick.label.set_fontsize(fontsize) for tick in ax2.yaxis.get_major_ticks()]

ax2.set_xlabel(r'growth rate $[h^{-1}]$', size=fontsize)
ax2.set_ylabel('median $r_{\mathrm{cat}}$ $[s^{-1}]$', size=fontsize)

plt.tight_layout()
plt.savefig('%s/pdf/rcat_increase_with_gr.pdf'%R.path)#