from kapp import *
import sys, os
#sys.path.append(os.path.expanduser('~/git/component-contribution'))
sys.path.append(os.path.expanduser('~/git/cobrapy'))
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)

rc = PLOT_DATA(model)

reactions = rc.calculate_enzyme_rates().index
gr = rc.growth_conditions.growth_rate_1_h
V = rc.V_data.loc[reactions]
E = rc.E_data.loc[reactions]

for c in combinations(gc, 2):
        
    i, j = sorted(c, key=gr.to_dict().get)

    flux_change = V[j]/V[i]
    expression_change = E[j]/E[i]
        
    relative_expression_change = expression_change / flux_change
    
    if (i, j) == ('acetate_heinmann', 'glucose_heinmann'):    
        x1 = relative_expression_change
        x1 = x1.replace([np.inf, -np.inf], np.nan).dropna()
    if (i, j) == ('Chemostat_vilu_011', 'Chemostat_vilu_049'):
        x2 = relative_expression_change
        x2 = x2.replace([np.inf, -np.inf], np.nan).dropna()

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(6,5), sharex=True, sharey=True)
ax1.hist(x1, histtype='stepfilled', bins=np.arange(-1, np.max(x1) + 0.1, 0.1), 
         color='#ff9c9f', edgecolor='none')
ax1.axvline(x1.median(), 0, 1, c=(0.8, 0, 0), zorder=10, ls='-', lw=3)
ax1.axvspan(0, 1, alpha=0.3, color='#cdb79e', zorder=3)

ax2.hist(x2, histtype='stepfilled', bins=np.arange(-1, np.max(x2) + 0.1, 0.1), 
         color='#ffe07e', edgecolor='none')
ax2.axvline(x2.median(), 0, 1, c=(0.6, 0.6, 0), zorder=10, ls='-', lw=3)
ax2.axvspan(0, 1, alpha=0.3, color='#cdb79e', zorder=3)

plt.xlim(-0.5,1.5)
plt.ylim(0,60)
plt.ylabel('reactions', size=15)

#ax1.set_title(r'acetate $(C_1)$ $\longrightarrow $ glucose $(C_2)$', weight='book', size=15)
#ax1.set_xlabel(r'fold change $\frac{E[{C_2}]}{E[{C_1}]}$', size=15)

#ax1.tick_params(axis='both', which='both', top='off', bottom='on',
#                                          left='on', right='off')
#ax1.set_ylim(0,60)
#ax1.grid(color='w', ls='-', lw=1, zorder=0)
#ax1.tick_params(color='w')

#plt.tight_layout()
#plt.savefig('../res/expression_regulation_acetate_to_glucose.svg')


#fig2 = plt.figure(figsize=(6,3))
#ax2 = plt.axes(axisbg=(0.95,0.92,0.90))
#rc.expression_plot(x2, y2, ax2)
#ax2.set_title(r'glucose limited chemosat - slow $(C_1)$ to fast $(C_2)$', weight='book', size=15)
#ax2.set_xlabel(r'fold change $\frac{E[{C_2}]}{E[{C_1}]}$', size=15)
#ax2.set_ylabel('reactions', size=15)   
#ax2.tick_params(axis='both', which='both', top='off', bottom='on',
#                                          left='on', right='off')
#ax2.set_ylim(0,60)
#ax2.grid(color='w', ls='-', lw=1, zorder=0)
#ax2.tick_params(color='w')
#        
#plt.tight_layout()#    


#ax5.tick_params(axis='y', which='both', left='on', right='off')
#ax6.tick_params(axis='y', which='both', left='off', right='off')
#plt.setp(ax6.get_yticklabels(), visible=False)


#index = []
#for i in growth_rates.index:
#    for j in growth_rates.index:
#        if growth_rates[j] / growth_rates[i] >= 1:
#            index.append('%s:%s'%(i,j))                
#
#data = pd.DataFrame(index=index, columns=['gr_change', 
#                                          'v_change', 
#                                          'e_change',
#                                          'e_explained'])                                   
#
#
#for c in index:
#    i, j = c.split(':')
#    
#    v2_to_v1 = np.log2(V[j] / V[i]).replace([np.inf, -np.inf], np.nan)
#    E2_to_E1 = np.log2(E[j] / E[i]).replace([np.inf, -np.inf], np.nan)
#    gr2_to_gr1 = np.log2(growth_rates[j] / growth_rates[i])
#
#    data['gr_change'][c] = gr2_to_gr1
#    data['v_change'][c] = v2_to_v1.dropna().median()                
#    data['e_change'][c] = E2_to_E1.dropna().median()
#    data['e_explained'][c] = (E2_to_E1 / v2_to_v1).dropna().median()
#
#    if c == 'acetate_heinmann:glucose_heinmann':    
#        x1, y1 = E2_to_E1,  v2_to_v1
#
#    if c == 'Chemostat_vilu_011:Chemostat_vilu_049':
#        x2, y2 = E2_to_E1, v2_to_v1 
#        
#fig1 = plt.figure(figsize=(6,3))
#ax1 = plt.axes(axisbg=(0.95,0.92,0.90))
#rc.expression_plot(x1, y1, ax1)
#ax1.set_title(r'acetate $(C_1)$ $\longrightarrow $ glucose $(C_2)$', weight='book', size=15)
#ax1.set_xlabel(r'fold change $\frac{E[{C_2}]}{E[{C_1}]}$', size=15)
#ax1.set_ylabel('reactions', size=15)   
#ax1.tick_params(axis='both', which='both', top='off', bottom='on',
#                                          left='on', right='off')
#ax1.set_ylim(0,60)
#ax1.grid(color='w', ls='-', lw=1, zorder=0)
#ax1.tick_params(color='w')
#
#plt.tight_layout()
#plt.savefig('../res/expression_regulation_acetate_to_glucose.svg')
#
#
#fig2 = plt.figure(figsize=(6,3))
#ax2 = plt.axes(axisbg=(0.95,0.92,0.90))
#rc.expression_plot(x2, y2, ax2)
#ax2.set_title(r'glucose limited chemosat - slow $(C_1)$ to fast $(C_2)$', weight='book', size=15)
#ax2.set_xlabel(r'fold change $\frac{E[{C_2}]}{E[{C_1}]}$', size=15)
#ax2.set_ylabel('reactions', size=15)   
#ax2.tick_params(axis='both', which='both', top='off', bottom='on',
#                                          left='on', right='off')
#ax2.set_ylim(0,60)
#ax2.grid(color='w', ls='-', lw=1, zorder=0)
#ax2.tick_params(color='w')
#        
#plt.tight_layout()
#plt.savefig('../res/expression_regulation_chemostat.svg')
