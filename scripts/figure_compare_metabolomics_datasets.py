import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import numpy as np

#metabolites = pd.DataFrame.from_csv("../data/model_metabolites_to_kegg.txt", sep='\t')
#metabolites = metabolites.kegg_id

palsson = pd.DataFrame.from_csv('../data/palsson_metabolite_concentrations[mM].csv')
bennett = pd.DataFrame.from_csv('../data/bennett_metabolite_concentrations[mM].csv')

cids = palsson.index & bennett.index

palsson = palsson['palsson_glc [mM]'][cids]
bennett = bennett['bennett_glc [mM]'][cids]

plt.figure(figsize=(6,6))
ax = plt.axes()
plt.scatter(palsson, bennett, alpha=0.5, edgecolor='none')
plt.plot([1e-3, 1e3], [1e-3, 1e3], 'r', zorder=0)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-3, 1e3)
ax.set_ylim(1e-3, 1e3)

ax.set_xlabel('McCloskey et al, [mM]', size=15)
ax.set_ylabel('Bennett et al, [mM]', size=15)

logx = np.log(palsson)
logy = np.log(bennett)
r, p = pearsonr(logx, logy)

ax.text(1e-3*5, 1e2, '$R^2$=%.2f' %r**2, size=15)

plt.tight_layout()
plt.savefig('../res/comparing_metabolomics.pdf')