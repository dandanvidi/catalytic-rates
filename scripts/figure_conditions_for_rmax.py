''' condition histogram ''' 
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
from rcat import RCAT
import sys, os
import numpy as np
sys.path.append(os.path.expanduser('~/git/across-projects'))
from color import ColorMap

R = RCAT()
rmax = R.rmax
conds = Counter(rmax.condition)

values = sorted(conds.values(), reverse=True)
labels = sorted(conds.keys(), key=conds.get, reverse=True)

fig = plt.figure(figsize=(6,6))
ax = plt.axes()

bar = plt.bar(range(len(conds)), values)
plt.xticks(np.arange(0,len(conds))+0.4, labels, rotation=90)

gc = pd.DataFrame.from_csv("../data/growth_conditions.csv")
cs = [gc['carbon source'][l] for l in labels]
cs_dist = Counter(cs)
colors = ColorMap(cs_dist)
for i, l in enumerate(labels):
    bar[i].set_color(colors[cs[i]])

plt.yticks(np.arange(0,41,10))
ax.set_ylabel('$r_\mathrm{max}$ values', size=15)
ax.set_ylim(0,41)
ax.tick_params(axis='both', which='both', top='off', right='off')
plt.tight_layout()

plt.savefig("../res/best_condition.svg")