from kapp import MM_KINETICS, PLOT_DATA
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)
mm = MM_KINETICS(model)
p = PLOT_DATA(model)

KMs = pd.DataFrame.from_csv('../cache/KM_values.csv')
kcat = mm.get_kcat_of_model_reactions()
rcat = mm.calculate_enzyme_rates()
rmax = mm.get_rcat_max(7)
cmax = mm.get_cond_max(7)
b_to_r = mm.map_model_reaction_to_genes()

S, T, ST = mm.concentraion_dependant_effects()

newrmax = rmax.div(ST).dropna()
reactions = kcat.index & rmax.index
corrected_reactions = newrmax.index & kcat.index
other_reactions = reactions - corrected_reactions

newrmax = newrmax[corrected_reactions]

fig = plt.figure(figsize=(6,6))
ax = plt.axes()

labels = []
for r in corrected_reactions:
    if newrmax[r] / rmax[r] > 1.5:
        labels.append(r)
    ax.vlines(kcat[r], rmax[r], newrmax[r], zorder=3, color='0.5')


ax.scatter(kcat[corrected_reactions],rmax[corrected_reactions], s=15, c='k', zorder=3, edgecolor='none')
report = p.plot_kcat_rcat_correlation(kcat[corrected_reactions], newrmax[corrected_reactions], fig, ax, 
                                color='#AA6939', yerr='none', labels=labels, fit_on=False)

                                
ax.set_ylabel(r'in vivo $r_{\mathrm{max}}$ / $ST$ $\left[s^{-1}\right]$', size=20, style='italic')
ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}$ $\left[s^{-1}\right]$', size=20)
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.ticklabel_format(size=40)

plt.tight_layout()
plt.savefig('../res/correaltion_after_ST.pdf')