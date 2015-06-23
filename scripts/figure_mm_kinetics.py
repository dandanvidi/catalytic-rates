from rcat import RCAT
from saturation_and_thermodynamics import MM_KINETICS
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from figure_correlation_plot import plot_kcat_rmaxn_correlation
from scipy import stats

R = RCAT()
mm = MM_KINETICS(R.model, list(R.rmax.index))
x = R.kcat['kcat [s^-1]'][mm.reactions]
y = R.rmaxn['rmax x n [s^-1]'][mm.reactions]

conditions = ['glc', 'ac', 'glyc']
s = pd.DataFrame(columns=conditions, index=mm.reactions)
t = pd.DataFrame(columns=conditions, index=mm.reactions)

for i, c in enumerate(conditions):
    s[c] = mm.get_saturation_effect(c)
    t[c] = mm.get_thermodynamic_effect(c)

invivo_kcat = pd.Series(index=mm.reactions)
invivo_kcat_uc = pd.Series(index=mm.reactions)
S = pd.Series(index=mm.reactions)
T = pd.Series(index=mm.reactions)
for r in y.index:
    c = R.rmaxn['carbon source'][r]
    rmx = R.rmaxn['rmax x n [s^-1]'][r]
    
    if c in conditions and s.loc[r][c]*t.loc[r][c]:
        S[r] = (s.loc[r][c]).n
        T[r] = (t.loc[r][c]).n
        invivo_kcat[r] = (rmx/(S[r]*T[r]))
        invivo_kcat_uc[r] = (rmx/(s.loc[r][c]*t.loc[r][c])).s
#        
reactions = invivo_kcat.dropna().index

y1 = invivo_kcat[reactions]
S = S[reactions]
T = T[reactions]
labels = {r:R.map_reactions_to_gene_names()[r] for r in reactions}
x = x[reactions]
#y_uc = kcat_vivo_uncertain[index]
#s = s[index]
#t = t[index]
#
fig = plt.figure(figsize=(6,6))
ax= plt.axes()
report1 = plot_kcat_rmaxn_correlation(x, y1, fig, ax, color='k', zorder=4)
rmse = np.sqrt( report1.sum_square / len(x) )
r1, pval1 = stats.pearsonr(np.log10(x), np.log10(y1))
##
##
##
def stacked_residual(x, y, saturation_effect, thermodynamic_effect, ax):
    for (x0, y0, s, t) in zip(x, y, saturation_effect, thermodynamic_effect):
        ax.vlines(x0, y0, y0/s, lw=3.5, colors='b', alpha=0.5)
        ax.vlines(x0, y0/s, y0/(s*t), lw=3.5, colors='#FFCC00', alpha=0.8)
##
y = y[reactions]
stacked_residual(x.values, y.values, S.values, T.values, ax)
##
report = plot_kcat_rmaxn_correlation(x, y, fig, ax, color='#AA6939', labels=labels, 
                                     hide_overlap=False)
rmse = np.sqrt( report.sum_square / len(x) )
r, pval = stats.pearsonr(np.log10(x), np.log10(y))
ax.set_ylabel(r'in vivo $r_{\mathrm{max}}\,\left[s^{-1}\right]$', 
              size=20, style='italic')
ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}\,\left[s^{-1}\right]$', 
              size=20, style='italic')
ax.tick_params(axis='both', which='both', top='off', right='off')

ax.set_xlim(1e-1*2,1e3*4)
ax.set_ylim(1e-1*2,1e3*4)
#
plt.tight_layout()
plt.savefig('%s/svg/saturation_and_thermodynamics_data.svg'%R.path)
#
