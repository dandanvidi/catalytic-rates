from rcat import RCAT
import pandas as pd
import matplotlib.pyplot as plt
from figure_correlation_plot import plot_kcat_rmaxn_correlation, add_labels

R = RCAT()

fva_0 = pd.DataFrame.from_csv('../cache/flux_variability[mmol_h_gCDW]_relaxation=0.csv')

fig = plt.figure(figsize=(6,6))
ax = plt.axes()

kcat = R.kcat['kcat [s^-1]']
rmaxn = R.rmaxn['rmax [s^-1]']
index = kcat.index & rmaxn.index
x = kcat[index]
y = rmaxn[index]

v = R.v.loc[index]
fba = pd.Series(index=index)
for r in index:
    c = R.rmaxn['condition'][r]
    fba[r] = v[c][r]

#    
#
#report = plot_kcat_rmaxn_correlation(x, y, 
#                                     fig, ax,
#                                     labels=labels, 
#                                     fit=True)
#
#
#plt.scatter(range(len(fva_0.index)), fva_0.minimum/fva_0.maximum, c='r', s=15)
#ax.set_xlim(0, len(fva_0.index))
#ax.set_ylim(0,1.5)
#
##max_rmaxn = 
#plot_kcat_rmaxn_correlation