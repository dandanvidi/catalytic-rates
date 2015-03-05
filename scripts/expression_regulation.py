from kapp import *
import sys, os
#sys.path.append(os.path.expanduser('~/git/component-contribution'))
sys.path.append(os.path.expanduser('~/git/cobrapy'))
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)

rc = PLOT_DATA(model)
growth_rates = rc.growth_conditions.growth_rate_1_h
    
reactions = rc.calculate_enzyme_rates().index

#ax5.tick_params(axis='y', which='both', left='on', right='off')
#ax6.tick_params(axis='y', which='both', left='off', right='off')
#plt.setp(ax6.get_yticklabels(), visible=False)


index = []
for i in growth_rates.index:
    for j in growth_rates.index:
        if growth_rates[j] / growth_rates[i] > 1:
            index.append('%s:%s'%(i,j))                

data = pd.DataFrame(index=index, columns=['gr_change', 
                                          'v_change', 
                                          'e_change',
                                          'e_explained'])                                   

V = rc.V_data.loc[reactions]
E = rc.E_data.loc[reactions]
for c in index:
    i, j = c.split(':')
    
    v2_to_v1 = np.log2(V[j] / V[i]).replace([np.inf, -np.inf], np.nan)
    E2_to_E1 = np.log2(E[j] / E[i]).replace([np.inf, -np.inf], np.nan)
    gr2_to_gr1 = np.log2(growth_rates[j] / growth_rates[i])

    data['gr_change'][c] = gr2_to_gr1
    data['v_change'][c] = v2_to_v1.dropna().median()                
    data['e_change'][c] = E2_to_E1.dropna().median()
    data['e_explained'][c] = (E2_to_E1 / v2_to_v1).dropna().median()

    if c == 'acetate_heinmann:glucose_heinmann':    
        x1, y1 = E2_to_E1,  v2_to_v1

    if c == 'Chemostat_vilu_011:Chemostat_vilu_049':
        x2, y2 = E2_to_E1, v2_to_v1 
        
fig1 = plt.figure(figsize=(6,3))
ax1 = plt.axes(axisbg='#FFE6C0')
rc.expression_plot(x1, y1, ax1)
ax1.set_title(r'acetate $\longrightarrow $ glucose', weight='book', color='b', size=15)
ax1.set_xlabel(r'fold change $\left[\frac{E[{C_1}]}{E[{C_2}]}\right]$', size=15)
ax1.set_ylabel('reactions', size=15)   
ax1.annotate(r'$\frac{\mu[{C_2}]}{\mu[{C_1}]} \approx \frac{v[{C_1}]}{v[{C_2}]}$', xy=(y1.mean(), 52),  xycoords='data',
            xytext=(y1.mean()+2.8, 48),
            horizontalalignment='right', verticalalignment='center',size=20
            )
ax1.tick_params(axis='both', which='both', top='off', bottom='on',
                                          left='on', right='off')
ax1.set_ylim(0,60)
plt.tight_layout()
plt.savefig('../res/expression_regulation_acetate_to_glucose.pdf')


fig2 = plt.figure(figsize=(6,3))
ax2 = plt.axes(axisbg='#FFE6C0')
rc.expression_plot(x2, y2, ax2)
ax2.set_title(r'glucose limited chemosat', weight='book', color='b', size=15)
ax2.set_xlabel(r'fold change $\left[\frac{E[{C_1}]}{E[{C_2}]}\right]$', size=15)
ax2.set_ylabel('reactions', size=15)   
ax2.annotate(r'$\frac{\mu[{C_2}]}{\mu[{C_1}]} \approx \frac{v[{C_2}]}{v[{C_1}]}$', xy=(y2.mean(), 52),  xycoords='data',
            xytext=(y2.mean()+0.05, 48),
            horizontalalignment='left', verticalalignment='center',size=20
            )
ax2.tick_params(axis='both', which='both', top='off', bottom='on',
                                          left='on', right='off')
ax2.set_ylim(0,60)
plt.tight_layout()
plt.savefig('../res/expression_regulation_chemostat.pdf')
