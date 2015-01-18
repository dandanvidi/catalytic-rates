from kapp import MODEL 
import scipy.odr
import math
import statsmodels.api as sm
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from cobra.io.sbml import create_cobra_model_from_sbml_file

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
M = MODEL(model)

def plot(x, y, fig, ax, color, edge, yerr=None, limits=None, labels=[]):

    logx = np.log10(x)
    logy = np.log10(y)
    
    if limits==None:
        a = math.ceil(logy.max())+2
        limits = np.array([10**(a-10), 10**a]) 
    
    # fit y = ax + b using robust linear method
#    rob_x = sm.add_constant(logx)
#    resrlm = sm.RLM(logy, rob_x).fit()
#    
#    intercept, slope = resrlm.params
    newx = logx
    newx = np.append(newx, 4.0)
    newx = np.append(-4.0, newx)
    
    ax.scatter(x,y,s=55, c=color, marker='o', edgecolor=edge)
    
    


#    ax.errorbar(x.PGI, y.PGI, 
#                yerr=yerr,
#                uplims=True,
#                barsabove=False,
#                fmt=None, 
#                ecolor='0.5', alpha=0.5)

    ax.plot([1e-4, 1e4], limits, 'k', ls='--')
#    ax.plot(10**newx, 10**(slope * np.array(newx) + intercept), color='#33CCFF',
#             label='Oridinary Least Squares')


    
    #Define function for scipy.odr
    fit_func = lambda B,x: B[0]*x + B[1]
    
    #Fit the data using scipy.odr
    Model = scipy.odr.Model(fit_func)
    Data = scipy.odr.RealData(np.log10(x), np.log10(y))
    Odr = scipy.odr.ODR(Data, Model, beta0=[1,1])
    output = Odr.run()
    #output.pprint()
    beta = output.beta
    betastd = output.sd_beta
    print beta, betastd
    ax.plot(10**newx, 10**fit_func(beta, newx), color='#FF0000')
    
    
#    ax.fill_between(10**newx, 
#                    10**(newx*(beta[0]+betastd[0])+(beta[1]-betastd[1])),
#                    10**(newx*(beta[0]-betastd[1])+(beta[1]+betastd[1])), 
#                    color='0.5', alpha=0.7)

    ax.set_xscale('log')
    ax.set_yscale('log', nonposy='clip')
    ax.legend(loc='upper left', fontsize=10)
    ax.set_xlim(limits[0:])
    ax.set_ylim(limits[0:])
    

    b_to_r = M.map_model_reaction_to_genes().set_index(0)
    
    if labels!=[]:
        ann = []
        for r in labels:
            ann.append(ax.text(x[r], y[r], M.gene_names[b_to_r.loc[r][1]]))
            
        mask = np.zeros(fig.canvas.get_width_height(), bool)
        
        fig.canvas.draw()
        
        for a in ann:
            bbox = a.get_window_extent()
            x0 = int(bbox.x0)
            x1 = int(math.ceil(bbox.x1))
            y0 = int(bbox.y0)
            y1 = int(math.ceil(bbox.y1))
        
            s = np.s_[x0:x1+1, y0:y1+1]
            if np.any(mask[s]):
                a.set_visible(False)
            else:
                mask[s] = True

    ax.grid()
    
    #    ax.set_xticks(np.logspace(-5,5,11), size=10)
    ax.set_xlim(1e-4,1e4)
    ax.set_ylim(limits[0],limits[1])
    #    plt.yticks(np.logspace(limits[0],limits[1],11), size=10)
    #    ax.set_yticklabels(10**np.linspace(-4,4,9), size=15)
    cor, pval = stats.pearsonr(logx, logy)
    
    ax.set_ylim(limits)
    rmse = np.sqrt( output.sum_square / len(x) )
    print "r^2 = %.3f, pval = %.2e"%(cor**2, pval)
    print "rmse = %.3f" % rmse    
    return output
