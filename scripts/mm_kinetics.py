from kapp import MM_KINETICS, PLOT_DATA
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)
mm = MM_KINETICS(model)

KMs = pd.DataFrame.from_csv('../cache/KM_values.csv')
kcat = mm.get_kcat_of_model_reactions()
rcat = mm.calculate_enzyme_rates()
rmax = mm.get_rcat_max(7)
cmax = mm.get_cond_max(7)

S, T, ST = mm.concentraion_dependant_effects()

temp = T.index & S.index

fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(8,3.5), sharey=True)
#ax1.set_axis_bgcolor((0.95,0.92,0.90))
#ax2.set_axis_bgcolor((0.95,0.92,0.90))
#ax3.set_axis_bgcolor((0.95,0.92,0.90))

ax1.scatter(S, (rmax/kcat)[temp], s=50, alpha=0.4, edgecolor='0.5', zorder=3)            
ax1.set_xlim(1/1000.0, 3)
ax1.set_ylim(1/1000.0, 30)
ax1.axhline(1, 1e-3,3, ls=':', lw=3, color='k')
ax1.axvline(1, 0, 1e3, ls=':', lw=3, color='k')
ax1.set_ylabel(r'$r_{\rm{max}}$ / $k_{\rm{cat}}$', size=20)
ax1.set_xlabel('$ S $', size=20)
ax1.tick_params(axis='both', which='both', top='off', right='off')
ax1.plot(np.logspace(-3,0), np.logspace(-3,0), 'r:', lw=3)

ax2.scatter(T, (rmax/kcat)[temp], s=50, color='#FFCC00', alpha=0.5, edgecolor='0.5', zorder=3)            
ax2.set_xlim(1/1000.0, 3)
ax2.set_ylim(1/1000.0, 30)
ax2.axhline(1, 1e-3,3, ls=':', lw=3, color='k')
ax2.axvline(1, 0, 1e3, ls=':', lw=3, color='k')
ax2.set_xlabel('$ T $', size=20)
ax2.tick_params(axis='both', which='both', top='off', right='off')
ax2.plot(np.logspace(-3,0), np.logspace(-3,0), 'r:', lw=3)

ax3.scatter(ST, (rmax/kcat)[ST.index], s=50, color='#00A352', alpha=0.4, edgecolor='0.5', zorder=3)            
ax3.set_xlim(1/1000.0, 3)
ax3.set_ylim(1/1000.0, 30)
ax3.axhline(1, 1e-3,3, ls=':', lw=3, color='k', zorder=0)
ax3.axvline(1, 0, 1e3, ls=':', lw=3, color='k', zorder=0)
ax3.set_xlabel(r'$ S $  x  $ T $', size=20)
ax3.tick_params(axis='both', which='both', top='off', right='off')
ax3.plot(np.logspace(-3,0), np.logspace(-3,0), 'r:', lw=3, zorder=0)

r_to_b = mm.map_model_reaction_to_genes()

names = []
for r in ST.index:
    names.append(mm.gene_names[r_to_b[r]])

labels = list(ST[ST<0.4].index)  + ['PGI', 'PFK', 'AMAOTr', 'PTPATi']

for l in labels:
    if ST[l]:
        ax1.scatter(S[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
        ax2.scatter(T[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
        ax3.scatter(ST[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
        name = mm.gene_names[r_to_b[l]]
        xy1 = (S[l], (rmax/kcat)[l])
        xy2 = (T[l], (rmax/kcat)[l])
        xy3 = (ST[l], (rmax/kcat)[l])
        xy4 = (1, (rmax/kcat)[l])
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
        
        if l not in ['AMAOTr', 'PTPATi']:
            ax1.scatter(S[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
            ax2.scatter(T[l], (rmax/kcat)[l], s=50, color='none',edgecolor='k', zorder=3)
            ax1.annotate(name, (S[l]/1.4, (rmax/kcat)[l]),
                     ha='right', va='center', size=12)
            ax2.annotate(name, (T[l]/1.4, (rmax/kcat)[l]),
                     ha='right', va='center', size=12)
        ax3.annotate(name, (ST[l]/1.4, (rmax/kcat)[l]),
                 ha='right', va='center', size=12)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax2.set_xscale('log')
ax2.set_yscale('log')

ax3.set_xscale('log')
ax3.set_yscale('log')

plt.tight_layout()

plt.savefig('../res/saturation_and_thermodynamics.pdf')

#
#


