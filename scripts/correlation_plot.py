''' catalytic rates '''

from kapp import RCAT, PLOT_DATA
import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)

R = RCAT(model)
growth_conditions = R.growth_conditions

kcat = R.get_kcat_of_model_reactions()
rcat = R.calculate_enzyme_rates()
rmax = R.get_rcat_max(7)

reactions = rmax.index & kcat.index

x = kcat[reactions]
y = rmax[reactions]

fig = plt.figure(figsize=(6,6))
ax = plt.axes(axisbg='0.95')


#xtpi = kcat.TPI
#ytpi = rcat.T.TPI.max()
#
#plt.plot(xtpi, ytpi, 'ro')

report = PLOT_DATA(model).plot_kcat_rcat_correlation(x, y, fig, ax, 
                                color='#FFB84D', yerr='none', labels=reactions)

ax.set_ylabel(r'in vivo $r^{\mathrm{max}}_{\mathrm{cat}}$ $\left[s^{-1}\right]$', size=20, style='italic')
ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}$ $\left[s^{-1}\right]$', size=20)
ax.tick_params(axis='both', which='both', top='off', right='off')

plt.tight_layout()
plt.savefig('../res/kcat_rmax_correlation.pdf')


''' find conditions with rmax values '''
fig = plt.figure(figsize=(8,5))
ax = plt.axes(axisbg='0.95')
ax2 = ax.twinx()

cmax = pd.DataFrame(rcat.idxmax(axis=1))
hist = cmax.apply(pd.value_counts)
bars = growth_conditions.join(hist)
bars.rename(columns={0:'cases'}, inplace=True)
bars.replace(np.nan, 0, inplace=True)
width = 0.8
ind = np.arange(0,len(bars)*2,2)


ax2.bar(ind + width/1.9, bars.cases, width=width, color='#FFB84D')
ax.bar(ind, bars.growth_rate_1_h, width=width, color='b', alpha=0.5)

ax.set_xticks(ind+width)
ax.set_xticklabels(bars.alternative_name, rotation='vertical')
ax.set_ylabel(r'growth rate $\left[h^{-1}\right]$', color='b', size=15)
ax2.set_ylabel('number of enzymes', color='#FFB84D', size=15)
ax.tick_params(axis='both', which='both', top='off', right='off', bottom='off')
ax.set_xlim([-0.5,len(bars)*2])
plt.tight_layout()
plt.savefig('../res/best_conditions.pdf')