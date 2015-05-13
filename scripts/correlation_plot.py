''' catalytic rates '''

from kapp import RCAT, PLOT_DATA
from cobra.io.sbml import create_cobra_model_from_sbml_file
import matplotlib.pyplot as plt

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)

R = RCAT(model)
growth_conditions = R.growth_conditions

b_to_r = R.map_model_reaction_to_genes()
kcat = R.get_kcat_of_model_reactions()
rcat = R.calculate_enzyme_rates()
rmax = R.get_rcat_max(7)

reactions = rmax.index & kcat.index

x = kcat[reactions]
y = rmax[reactions]

fig = plt.figure(figsize=(6,6))
ax = plt.axes()

if 'TPI' in x.index:
    ax.text(x.TPI*1.2, y.TPI, R.gene_names[b_to_r['TPI']], 
                           ha='left', va='center')

report = PLOT_DATA(model).plot_kcat_rcat_correlation(x, y, fig, ax, 
                                color='#AA6939', yerr='none', labels=list(reactions), fit_on=False)

ax.set_ylabel(r'in vivo $r_{\mathrm{max}}$ $\left[s^{-1}\right]$', size=20, style='italic')
ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}$ $\left[s^{-1}\right]$', size=20)
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.ticklabel_format(size=40)

plt.tight_layout()
plt.savefig('../res/kcat_rmax_correlation.svg')
