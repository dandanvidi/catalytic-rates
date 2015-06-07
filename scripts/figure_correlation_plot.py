from rcat import RCAT
from scipy import stats, odr
import numpy as np
import matplotlib.pyplot as plt
import math

def plot_kcat_rcat_correlation(x, y, fig, ax, color='b', edge='none', 
                               yerr='none', labels={}, hide_overlap=True,
                               fit=False, zorder=3):
    
    logx = np.log10(x)
    logy = np.log10(y)
    
    ax.scatter(x, y,s=20, c=color, alpha=0.8, marker='o', edgecolor=edge, zorder=zorder)
    
    if yerr != 'none':
        ax.errorbar(x, y, 
                    yerr=yerr, barsabove=False, 
                    fmt=None, ecolor='k', alpha=0.4)
                
    ax.plot([1e-4, 1e4], [1e-4,1e4], '#333676', ls='-', lw=2, zorder=5)
     
    #Define function for scipy.odr
    fit_func = lambda B,x: B[0]*x + B[1]
    
    #Fit the data using scipy.odr
    Model = odr.Model(fit_func)
    Data = odr.RealData(logx, logy)
    Odr = odr.ODR(Data, Model, beta0=[1,1])
    output = Odr.run()
    #output.pprint()
    beta = output.beta                

    if fit:  
        edge = np.array([-4, 4])
        ax.plot([1e-4, 1e4], 10**fit_func(beta, edge), color='#699A33', ls=':', lw=3, zorder=1)
        
        
    ax.set_xscale('log', nonposx='clip')
    ax.set_yscale('log', nonposy='clip')
                    
    if labels!={}:
        ann = []
        for r, name in labels.iteritems():
            if x[r]>y[r]:
                ann.append(ax.text(x[r], y[r]/1.1, name, 
                                    ha='center', va='top', zorder=5))
            if x[r]<y[r]:
                ann.append(ax.text(x[r], y[r]*1.1, name,
                                    ha='center', va='bottom', zorder=5))
                                        
            mask = np.zeros(fig.canvas.get_width_height(), bool)
            fig.canvas.draw()
            for i, a in enumerate(ann):
                bbox = a.get_window_extent()
                x0 = int(bbox.x0)
                x1 = int(math.ceil(bbox.x1))
                y0 = int(bbox.y0)
                y1 = int(math.ceil(bbox.y1))
            
                s = np.s_[x0:x1, y0:y1]
                if hide_overlap:
                    if np.any(mask[s]):
                        a.set_visible(False)
                    else:
                        mask[s] = True
                else:
                    mask[s] = True
                
                
    for tickx, ticky in zip(ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()):
        tickx.label.set_fontsize(14) 
        ticky.label.set_fontsize(14)
        
    ax.set_xlim(1e-3/4,4*1e3)
    ax.set_ylim(1e-3/4,4*1e3)

    return output

if __name__ == "__main__":
    
    R = RCAT()
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes()
    rcat = R.rcat
    kcat = R.kcat
    rmax = R.rmax['rmax [s^-1]']
    index = kcat.index & rmax.index
    x = kcat[index]
    y = rmax[index]
    labels = {k:v for k,v in R.map_reactions_to_gene_names().iteritems() if k in index}
    report = plot_kcat_rcat_correlation(x, y, fig, ax, 
                               color='#AA6939', edge='none', 
                               yerr='none', labels=labels, 
                               fit=False)
    
    
    rmse = np.sqrt( report.sum_square / len(x) )
    r, pval = stats.pearsonr(np.log10(x), np.log10(y))
    ax.text(1e-3/2, 1e2/1.25, '$R^2=$%.2f\n$RMSE=$%.2f' %(r**2,rmse), size=15)
    
    # specif labels to add
    ax.text(x['TPI'], y['TPI']*1.1, 'tpi', ha='center', va='bottom', zorder=5)
    ax.text(x['PFK'], y['PFK']*1.1, 'pfkA', ha='center', va='bottom', zorder=5)
    
    ax.set_ylabel(r'in vivo $r_{\mathrm{max}}$ $\left[s^{-1}\right]$', 
                  size=20, style='italic')
    ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}$ $\left[s^{-1}\right]$', 
                  size=20, style='italic')
    ax.tick_params(axis='both', which='both', top='off', right='off')
    ax.ticklabel_format(size=40)

    ax.set_xlim(1e-3/4,4*1e3)
    ax.set_ylim(1e-3/4,4*1e3)
    
    plt.tight_layout()
    plt.savefig('../res/kcat_rmax_correlation.svg')
