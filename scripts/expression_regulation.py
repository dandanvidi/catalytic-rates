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
from matplotlib import cm

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
rc = RCAT(model)

reactions = rc.calculate_enzyme_rates().index
gc = rc.growth_conditions
gm = gc.growth_mode
V = rc.V_data.loc[reactions]
V.replace([0.0, -0.0], np.nan, inplace=True)
E = rc.E_data.loc[reactions]

cmap = cm.get_cmap('Reds')

def draw_fold_change_hist(ax, array, color):
    ax.hist(array, histtype='stepfilled', 
            bins=np.arange(-1, np.max(array) + 0.1, 0.1),                                                         
            color=color, edgecolor='none')
    ax.axvline(array.median(), 0, 1, c='k', zorder=10, ls='-', lw=3)
    ax.axvspan(array.median(), 1, alpha=0.3, color='#cdb79e', zorder=3)
    ax.axvline(1, 0, 1, c='k', zorder=10, ls=':', lw=3)
    
def get_flux_and_expression_fold_change(growthconditions, expression, fluxes):
    growthrates = growthconditions.growth_rate_1_h
    reactions = expression.index & fluxes.index
    combs = [(i,j) for (i,j) in combinations(growthrates.index, 2) 
                                    if growthrates[j] > growthrates[i]]
    mu_ratio = np.array([growthrates[j]/growthrates[i] for (i,j) in combs])
    mu2_to_mu1 = pd.DataFrame(index=combs, columns=['mu_ratio'])
    mu2_to_mu1.mu_ratio = mu_ratio

    E2_to_E1 = pd.DataFrame(index=reactions, columns=combs)
    V2_to_V1 = pd.DataFrame(index=reactions, columns=combs)    
    for i,j in combs:

        V2_to_V1[(i,j)] = np.log(fluxes[j]/fluxes[i])
        E2_to_E1[(i,j)] = np.log(expression[j]/expression[i])
        
    return E2_to_E1, V2_to_V1, mu2_to_mu1.mu_ratio


GC = gc
E2_to_E1, V2_to_V1, mu2_to_mu1 = get_flux_and_expression_fold_change(GC, E, V)

E2E1_to_V2V1 = E2_to_E1/V2_to_V1
E2E1_to_V2V1.replace([np.inf, -np.inf], np.nan, inplace=True)
E2E1_to_V2V1.dropna(how='any', inplace=True)

x = mu2_to_mu1
y = E2E1_to_V2V1.median()
print y[('Chemostat_vilu_011', 'Chemostat_vilu_049')]
print E2_to_E1[('Chemostat_vilu_011', 'Chemostat_vilu_049')].median()
print V2_to_V1[('Chemostat_vilu_011', 'Chemostat_vilu_049')].median()
weighted_by_gr = y*(mu2_to_mu1-1)
explained = weighted_by_gr.sum()/(mu2_to_mu1-1).sum()
print explained

#GC = gc[gc.growth_mode=='chemostat']
#E2_to_E1, V2_to_V1, mu2_to_mu1 = get_flux_and_expression_fold_change(GC, E, V)
#
#E2E1_to_V2V1 = E2_to_E1/V2_to_V1
#E2E1_to_V2V1.replace([np.inf, -np.inf], np.nan, inplace=True)
#E2E1_to_V2V1.dropna(how='any', inplace=True)
#
x = mu2_to_mu1
#y = E2E1_to_V2V1.median()*(mu2_to_mu1-1)
#plt.scatter(x,y, alpha=0.5)

plt.scatter(V2_to_V1.median(),E2_to_E1.median(), alpha=0.5)
plt.plot(V2_to_V1.median(),V2_to_V1.median())
#plt.hlines(y.median(), 0, 7)
#plt.xlim(0,7)
#plt.ylim(-1,2)

#weighted_by_gr = E2E1_to_V2V1.median()*(mu2_to_mu1-1)
#explained = weighted_by_gr.sum()/mu2_to_mu1.sum()
#print explained

#plt.figure()
#for (i,j) in E2E1_to_V2V1.columns:
#    y = E2E1_to_V2V1[(i,j)]
#    plt.scatter(range(len(y)), y, alpha=0.5)
#    print y.median()*mu2_to_mu1[(i,j)]


#plt.figure()
#plt.scatter(E2E1_to_V2V1)    
#i,j = ('acetate_heinmann', 'glucose_heinmann')
#y = E2E1_to_V2V1[(i,j)]
#y.replace([np.inf, -np.inf], np.nan, inplace=True)
#y.dropna(inplace=True)
#plt.scatter(range(len(y)), y, alpha=0.5, c='r')
#print y.mean(), y.median()1
#i,j = ('pyruvate_heinmann', 'glycerol_heinmann')
#y = (E2_to_E1[(i,j)]/V2_to_V1[(i,j)])
#y.replace([np.inf, -np.inf], np.nan, inplace=True)
#y.dropna(inplace=True)
#plt.scatter(range(len(y)), y, alpha=0.5, c='y')

#plt.xscale('log')
#plt.yscale('log')
#plt.scatter(range(len(y2)), y2, color='r')

#for i, r in enumerate(reactions):
#    if (y1.median()+0.01) <= y[r] or y[r] <= (y1.median()-0.01):
#        print y[r]
#        plt.text(i, y1[r], rc.gene_names[r_to_b[r]])
#
#fig, (ax1, ax2) = plt.subplots(2,1, figsize=(4,5), sharex=True, sharey=True)
#
#pairwise_E_to_V = pd.DataFrame(index=['median_fold_change','gr_ratio'])
#E_to_V = pd.DataFrame(index=gc.index, columns=gc.index)
#for c in combinations_with_replacement(gc.index, 2):        
#    i, j = sorted(c, key=gr.to_dict().get)
#    if i == j:
#        flux_change = np.nan
#    else:
#        flux_change = V[j]/V[i]
#    expression_change = E[j]/E[i]        
#    relative_expression_change = expression_change / flux_change
#    x = relative_expression_change
#    x = x.replace([np.inf, -np.inf], np.nan).dropna()
#    E_to_V[j][i] = x.median()
#    pairwise_E_to_V[c] = 1.
#    pairwise_E_to_V[c]['median_fold_change'] = x.median()
#    pairwise_E_to_V[c]['gr_ratio'] = gr[j]/gr[i]
#    if (i,j) == ('acetate_heinmann','glucose_heinmann'):
#        draw_fold_change_hist(ax1, x, cmap(x.median()))
#        print x.median()
#    elif (i,j) == ('Chemostat_vilu_011','Chemostat_vilu_049'):    
#        draw_fold_change_hist(ax2, x, cmap(x.median()))
#        print x.median()
#
#plt.xlim(0,1.5)
#plt.ylim(0,60)
#plt.ylabel('reactions', size=15)
#
#plt.savefig('../res/expression_regulation_hists.svg')
#
#plt.figure(figsize=(5,5))
#E_to_V.replace(np.nan, -np.inf, inplace=True)
#E_to_V[E_to_V>1] = 1
#ax = plt.axes()
#
#cmap.set_under('w')
#im = imshow(E_to_V, interpolation="nearest", cmap=cmap)
#ax.set_xticks([])
#ax.set_yticks([])
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.colorbar(im, cax=cax)
#plt.clim(0,1)
#plt.tight_layout()
#plt.savefig('../res/expression_regulation_matrix.svg')
#
#plt.figure()
#pairwise_E_to_V = pairwise_E_to_V.T.sort(columns='gr_ratio').dropna(how='any')
#x = pairwise_E_to_V.gr_ratio
#y = pairwise_E_to_V.median_fold_change
#plt.scatter(x,y)
#plt.xlabel('$\mu(C_H) / \mu(C_L)$')
#plt.ylabel('$\E(C_H) / \mu(C_L)$')
#plt.savefig('../res/expression_regulation_by_growthrate.svg')
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
