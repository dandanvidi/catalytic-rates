from kapp import RCAT
from numpy import ceil
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from pandas.stats.api import ols
import numpy as np
import sys, os
import matplotlib.pyplot as plt
sys.path.append(os.path.expanduser('~/git/global_functions'))
from color import ColorMap
from matplotlib import rc
from scipy.stats import pearsonr
# activate latex text rendering
rc('text', usetex=True)

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)

R = RCAT(model)
r_to_b = R.map_model_reaction_to_genes()

GR = R.growth_conditions.growth_rate_1_h
rcat = R.calculate_enzyme_rates()
rmax = R.get_rcat_max(7)
kcat = R.get_kcat_of_model_reactions()
reactions = rmax.index & kcat.index
rcat = rcat.loc[reactions]
names = {r:R.gene_names[r_to_b[r]] for r in reactions}
reactions = sorted(names.keys(), key=names.get)
colors = ColorMap(reactions)

fig, axarr = plt.subplots(int(ceil(len(reactions)/4.0)),4, 
                          figsize=(15,60), sharex=True, sharey=True)
j = 0 
for i, r in enumerate(reactions):
    y = rcat.loc[r].dropna()
    y = y / y.mean()
    x = GR[y.index]
    
    rsq = pearsonr(x,y)[0]**2
    results = ols(y=y, x=x)
    a = results.beta.x
    b = results.beta.intercept
    g = names[r]
    axarr[j/4, j%4].set_axis_bgcolor('#FFE6C0')
    axarr[j/4, j%4].plot(x, y, c=colors[r], marker='o', lw=0, zorder=3)
    axarr[j/4, j%4].set_ylim(0,2.5)  
    axarr[j/4, j%4].set_xlim(0,0.7)
    
    axarr[j/4, j%4].grid(color='w', ls='-', lw=1, zorder=0)
    axarr[j/4, j%4].tick_params(color='w')        
    axarr[j/4, j%4].set_title(r"\textit{\textbf{%s}} [ y=%.1fx%+.1f; $r^2$=%.2f ]" %(g, a, b, rsq))

    j += 1

for i in range(4):
    axarr[(j-1)/4, i].set_xlabel('growth rate [$h^{-1}$]')
for i in range((j-1)/4):
    axarr[i, 0].set_ylabel('$r_{\mathrm{cat}}$ normalized')
#plt.figtext(0.4,0.06,'growth rate $h^{-1}$', size=35)
#plt.figtext(0.04,0.7,'relative change in $r_\rm{cat}$ $s^{-1}$', size=35, rotation='vertical')
plt.tight_layout()
plt.savefig('../res/single_enzyme_rate_change.pdf')

