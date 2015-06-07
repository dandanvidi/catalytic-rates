from rcat import RCAT
from saturation_and_thermodynamics import MM_KINETICS
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from figure_correlation_plot import plot_kcat_rcat_correlation
from scipy import stats

R = RCAT()
kcat = R.kcat
rcat = R.rcat
mm = MM_KINETICS(R.model, list(R.rmax.index))
kcat = kcat[mm.reactions]
rmax = R.rmax.loc[mm.reactions]
rmax['carbon source'] = map(lambda x: R.gc['carbon source'][x], rmax.condition)

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
for r in rmax.index:
    c = rmax['carbon source'][r]
    rmx = rmax['rmax [s^-1]'][r]
    
    if c in conditions and s.loc[r][c]*t.loc[r][c]:
        S[r] = (s.loc[r][c]).n
        T[r] = (t.loc[r][c]).n
        invivo_kcat[r] = (rmx/(S[r]*T[r]))
        invivo_kcat_uc[r] = (rmx/(s.loc[r][c]*t.loc[r][c])).s
        
reactions = invivo_kcat.dropna().index

y = invivo_kcat[reactions]
S = S[reactions]
T = T[reactions]
labels = {r:R.map_reactions_to_gene_names()[r] for r in reactions}
x = kcat[reactions]
#y_uc = kcat_vivo_uncertain[index]
#s = s[index]
#t = t[index]
#
fig = plt.figure(figsize=(6,6))
ax= plt.axes()
report = plot_kcat_rcat_correlation(x, y, fig, ax, color='k', zorder=4)
rmse = np.sqrt( report.sum_square / len(x) )
r, pval = stats.pearsonr(np.log10(x), np.log10(y))
#
#
#
def stacked_residual(x, y, saturation_effect, thermodynamic_effect, ax):
    for (x0, y0, s, t) in zip(x, y, saturation_effect, thermodynamic_effect):
        ax.vlines(x0, y0, y0/s, lw=3.5, colors='b', alpha=0.5)
        ax.vlines(x0, y0/s, y0/(s*t), lw=3.5, colors='#FFCC00', alpha=0.8)
#
#
y1 = rmax['rmax [s^-1]'][reactions]
stacked_residual(x.values, y1.values, S.values, T.values, ax)
#
report1 = plot_kcat_rcat_correlation(x, y1, fig, ax, color='#AA6939', labels=labels, 
                                     hide_overlap=False)
rmse1 = np.sqrt( report1.sum_square / len(x) )
r1, pval1 = stats.pearsonr(np.log10(x), np.log10(y1))
#
ax.text(1e-0/5, 1e2*5, '$R^2=$%.2f (%.2f)' %(r**2,r1**2), size=15)
ax.set_ylabel(r'in vivo $r_{\mathrm{max}}$ $\left[s^{-1}\right]$', 
              size=20, style='italic')
ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}$ $\left[s^{-1}\right]$', 
              size=20, style='italic')
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.ticklabel_format(size=40)

ax.set_xlim(1e-1,1e3*2)
ax.set_ylim(1e-1,1e3*2)
#
plt.tight_layout()
plt.savefig('../res/saturation_and_thermodynamics.svg')

#cmax = mm.get_cond_max(7)
#
#S, T, ST = mm.concentraion_dependant_effects()
#
#temp = T.index & S.index
#newtemp = ((rmax/ST)[temp].dropna()).index
#x = kcat[newtemp]
#y = rmax[newtemp]
#S = S[newtemp]
#T = T[newtemp]
#ST = ST[newtemp]
#
#fig = plt.figure(figsize=(6,6))
#ax = plt.axes()
#report = pld.plot_kcat_rcat_correlation(x, y, fig, ax, color='#AA6939', edge='none', 
#                                   yerr='none', labels=[], fit_on=False)
#
#report = pld.plot_kcat_rcat_correlation(x, (y/ST)[newtemp], fig, ax, color='#333333', edge='none', 
#                                   yerr='none', labels=list(newtemp), fit_on=False)
#                                   
#def stacked_residual(x, y, S, T, ax):
##    resid = [sorted([y0,x0]) for (y0,x0) in zip(y,x)]
#    decom = [(1/s, 1/s/t) for (s, t) in zip(S, T)]
#    for x0, y0, (s,t) in zip(x, y, decom):
#        ax.vlines(x0, y0, y0*s, lw=3.5, colors='b', alpha=0.5)
#        ax.vlines(x0, y0*s, y0*t, lw=3.5, colors='#FFCC00', alpha=0.5)
##        if y0*t < x0:
##            ax.vlines(x0, y0*t, x0, lw=3.5, colors='#00B85C', linestyles='dotted', alpha=0.5)
##        if y0 > x0:
##            ax.vlines(x0, x0, y0, lw=3.5, colors='#00B85C', linestyles='dotted', alpha=0.5)
#
##        else:
##            ax.vlines(x0, x0, t, lw=3, colors='r', alpha=0.5)
#stacked_residual(x.values, y.values, S.values, T.values, ax)
#
#cor, pval = stats.pearsonr(np.log(x), np.log(y))
#rsq = cor**2
#newy = (y/ST).dropna()
#
#cor, pval = stats.pearsonr(np.log(x[newy.index]), np.log(newy))
#newrsq = cor**2
#ax.text(1e-2/2, 1e3, '$R^2=$%.2f (%.2f)' %(newrsq, rsq), size=15) 
#
#
##plt.rc('font', family='serif')
#ax.set_ylabel(r'in vivo $r_{\mathrm{max}}\,\left[s^{-1}\right]$', size=20, style='italic')
#ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}\,\left[s^{-1}\right]$', size=20, style='italic')
#ax.set_xlim(1e-3, 1e4)
#ax.set_ylim(1e-3, 1e4)
#ax.tick_params(axis='both', which='both', top='off', right='off')
#ax.ticklabel_format(size=40)
#
#plt.tight_layout()
#plt.savefig('../res/stacked_resid.pdf')
#fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(8,3.5), sharey=True)
##ax1.set_axis_bgcolor((0.95,0.92,0.90))
##ax2.set_axis_bgcolor((0.95,0.92,0.90))
##ax3.set_axis_bgcolor((0.95,0.92,0.90))
#
#ax1.scatter(S, (rmax/kcat)[temp], s=50, alpha=0.4, edgecolor='0.5', zorder=3)            
#ax1.set_xlim(1/1000.0, 3)
#ax1.set_ylim(1/1000.0, 30)
#ax1.axhline(1, 1e-3,3, ls=':', lw=3, color='k')
#ax1.axvline(1, 0, 1e3, ls=':', lw=3, color='k')
#ax1.set_ylabel(r'$r_{\rm{max}}$ / $k_{\rm{cat}}$', size=20)
#ax1.set_xlabel('$ S $', size=20)
#ax1.tick_params(axis='both', which='both', top='off', right='off')
#ax1.plot(np.logspace(-3,0), np.logspace(-3,0), 'r:', lw=3)
#
#ax2.scatter(T, (rmax/kcat)[temp], s=50, color='#FFCC00', alpha=0.5, edgecolor='0.5', zorder=3)            
#ax2.set_xlim(1/1000.0, 3)
#ax2.set_ylim(1/1000.0, 30)
#ax2.axhline(1, 1e-3,3, ls=':', lw=3, color='k')
#ax2.axvline(1, 0, 1e3, ls=':', lw=3, color='k')
#ax2.set_xlabel('$ T $', size=20)
#ax2.tick_params(axis='both', which='both', top='off', right='off')
#ax2.plot(np.logspace(-3,0), np.logspace(-3,0), 'r:', lw=3)
#
#ax3.scatter(ST, (rmax/kcat)[ST.index], s=50, color='#00A352', alpha=0.4, edgecolor='0.5', zorder=3)            
#ax3.set_xlim(1/1000.0, 3)
#ax3.set_ylim(1/1000.0, 30)
#ax3.axhline(1, 1e-3,3, ls=':', lw=3, color='k', zorder=0)
#ax3.axvline(1, 0, 1e3, ls=':', lw=3, color='k', zorder=0)
#ax3.set_xlabel(r'$ S $  x  $ T $', size=20)
#ax3.tick_params(axis='both', which='both', top='off', right='off')
#ax3.plot(np.logspace(-3,0), np.logspace(-3,0), 'r:', lw=3, zorder=0)
#
#r_to_b = mm.map_model_reaction_to_genes()
#
#names = []
#for r in ST.index:
#    names.append(mm.gene_names[r_to_b[r]])
#
#labels = list(ST[ST<0.4].index)  + ['PGI', 'PFK', 'AMAOTr', 'PTPATi']
#
#for l in labels:
#    if ST[l]:
#        ax1.scatter(S[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
#        ax2.scatter(T[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
#        ax3.scatter(ST[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
#        name = mm.gene_names[r_to_b[l]]
#        xy1 = (S[l], (rmax/kcat)[l])
#        xy2 = (T[l], (rmax/kcat)[l])
#        xy3 = (ST[l], (rmax/kcat)[l])
#        xy4 = (1, (rmax/kcat)[l])
#        con = ConnectionPatch(xyA=xy2, xyB=xy1, 
#                              coordsA="data", coordsB="data", 
#                              axesA=ax2, axesB=ax1, 
#                              color='0.5', alpha=0.5,
#                              shrinkA=4.25, shrinkB=4.25, lw=0.5)
#        ax2.add_artist(con)
#        
#        con = ConnectionPatch(xyA=xy3, xyB=xy2, 
#                      coordsA="data", coordsB="data", 
#                      axesA=ax3, axesB=ax2, 
#                      color='0.5', alpha=0.5,
#                      shrinkA=4.25, shrinkB=4.25, lw=0.5)
#        ax3.add_artist(con)

#        ax3.plot((ST[l]*1.25, 10), ((rmax/kcat)[l], (rmax/kcat)[l]), 
#                 c='0.5', alpha=0.5, lw=0.5)
        
#        if l not in ['AMAOTr', 'PTPATi']:
#            ax1.scatter(S[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
#            ax2.scatter(T[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
#            ax1.annotate(name, (S[l]/1.4, (rmax/kcat)[l]),
#                     ha='right', va='center', size=12)
#            ax2.annotate(name, (T[l]/1.4, (rmax/kcat)[l]),
#                     ha='right', va='center', size=12)
#        ax3.annotate(name, (ST[l]/1.4, (rmax/kcat)[l]),
#                 ha='right', va='center', size=12)
#
#ax1.set_xscale('log')
#ax1.set_yscale('log')
#
#ax2.set_xscale('log')
#ax2.set_yscale('log')
#
#ax3.set_xscale('log')
#ax3.set_yscale('log')
#
#plt.tight_layout()
#
#plt.savefig('../res/saturation_and_thermodynamics.pdf')
#
#res = kcat/rmax
#vitro = res[res>1]
#vivo = res[res<1]
#
#


