from rcat import RCAT
import numpy as np
import matplotlib.pyplot as plt

def boot_strap(x):
    x.dropna(inplace=True)
    median_values = np.zeros(1000)
    for i in xrange(1000):
        x = x.dropna()
        new_x = np.random.choice(x, len(x), replace=True)
        median_values[i] = np.median(new_x)
    return np.std(median_values)


R = RCAT()

rcat = R.rcat
kcat = R.kcat
rmax = R.rmax

reactions = rmax.index & kcat.index

rcat = rcat.loc[reactions]
rmax = rmax.loc[reactions]
kcat = kcat.loc[reactions]

gr = R.gc['growth rate (h-1)']
med = rcat.dropna(how='all').median()[R.gc.index]

grchemo = gr[R.gc[R.gc['growth mode']=='chemostat'].index]
medchemo = med[R.gc[R.gc['growth mode']=='chemostat'].index]


grbatch = gr[R.gc[R.gc['growth mode']=='batch'].index]
medbatch = med[R.gc[R.gc['growth mode']=='batch'].index]

fig = plt.figure(figsize=(6,5))
ax = plt.axes()

ax.scatter(grchemo, medchemo, zorder=10, c='#324EAC', 
           edgecolor='none', s=50, alpha=0.75, label='chemostat')
ax.scatter(grbatch, medbatch, zorder=10, c='#3FBE57', 
           edgecolor='none', s=50, alpha=0.75, label='batch')

ax.axhline(kcat.median(), color='k', ls=':', lw=3)
for i, c in enumerate(gr.index):
    if c in grbatch.index:
        ax.errorbar(gr[c], med[c], boot_strap(rcat[c]), c='#3FBE57', zorder=10, alpha=0.5)
    else:        
        ax.errorbar(gr[c], med[c], boot_strap(rcat[c]), c='#324EAC', zorder=10, alpha=0.5)

ax.text(.1, kcat.median()-0.25, 
         'median $k_{\mathrm{cat}}$', 
         va='top', size=15)
ax.set_xticks(np.arange(0,0.8,0.2))
ax.set_yticks(np.arange(0,10.1,2))
ax.set_xticklabels(np.arange(0,0.8,0.2), size=12.5)
ax.set_yticklabels(np.arange(0,11,2), size=12.5)
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.set_xlabel(r'growth rate $[h^{-1}]$', size=20)
ax.set_ylabel('median $r_{\mathrm{cat}}$ $[s^{-1}]$', size=20)
plt.legend(loc=1, fontsize=15, scatterpoints=1)
ax.ticklabel_format(size=40)
plt.tight_layout()
plt.savefig('../res/rcat_increase_with_gr.svg')#