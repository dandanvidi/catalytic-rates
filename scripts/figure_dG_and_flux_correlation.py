'''Radhakrishnan Mahadevan'''

from rcat import RCAT
from saturation_and_thermodynamics import MM_KINETICS
import matplotlib.pyplot as plt
R = RCAT()

mm = MM_KINETICS(R.model, list(R.rcat.index))

dGm = dict(zip(mm.reactions, mm.get_udGm_prime()))

rcat = R.rcat['Glucose'][mm.reactions].dropna()
dGm = {k:v for k, v in dGm.iteritems() if k in rcat.index}


x = [dGm[r].n for r in rcat.index]
y = [rcat[r] for r in rcat.index]

plt.figure()
ax = plt.axes()
plt.scatter(x,y, c='b', edgecolor='none', alpha=0.6)

ax.set_xlim(-200, 50)
ax.set_ylim(1e-4, 1e2*7)

ax.set_xlabel('$\Delta G\'_m \, [kJ/mol]$', size=15)
ax.set_ylabel('flux per enzyme $[s^{-1}]$', size=15)
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.set_yscale('log')

plt.tight_layout()
plt.savefig('../res/correlation_between_dG_and_flux.pdf')