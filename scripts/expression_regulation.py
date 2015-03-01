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
        
fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, sharey=True, figsize=(4,5))
#ax1 = fig.add_subplot(211, axis_bgcolor=(0.95,0.92,0.90))
#ax2 = fig.add_subplot(212, axis_bgcolor=(0.95,0.92,0.90))

rc.expression_plot(x1, y1, ax1)
rc.expression_plot(x2, y2, ax2)


ax1.set_title(r'acetate $\longrightarrow $ glucose', weight='book', color='b')
ax2.set_title(r'glucose limited chemostats', weight='book', color='b')

ax1.set_ylabel('reactions', size=15)   
ax2.set_ylabel('reactions', size=15)   

ax1.tick_params(axis='both', which='both', top='off', bottom='on',
                                          left='on', right='off')

ax2.tick_params(axis='both', which='both', top='off', bottom='on',
                                          left='off', right='off')
#ax2.set_yticks([])
ax2.set_xlabel(r'fold change $\left[\frac{E_n}{E_m}\right]$', size=15)
plt.tight_layout()

plt.savefig('../res/expression_regulation.pdf')
