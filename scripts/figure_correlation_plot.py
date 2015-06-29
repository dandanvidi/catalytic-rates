from rcat import RCAT
from scipy import stats, odr
import numpy as np
import matplotlib.pyplot as plt
import math

def plot_kcat_rmaxn_correlation(x, y, fig, ax, color='#9E7E5E', edge='none', 
                               yerr='none', labels={}, scatter_size=30, hide_overlap=True,
                               fit=False, zorder=3):
    
    logx = np.log10(x)
    logy = np.log10(y)
    
    ax.scatter(x, y,s=scatter_size, c=color, marker='o', edgecolor=edge, zorder=zorder)
    
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
        add_labels(x, y, labels, ax, fig, hide_overlap)
                                
    return output

def add_labels(x, y, labels, ax, fig, hide_overlap=True):
    ann = []
    for r, name in labels.iteritems():
        if x[r]>y[r]:
            ann.append(ax.text(x[r], y[r]/1.1, name, 
                                ha='center', va='top', zorder=5, size=13))
        if x[r]<y[r]:
            ann.append(ax.text(x[r], y[r]*1.1, name,
                                ha='center', va='bottom', zorder=5, size=13))
                                    
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

if __name__ == "__main__":
    
    R = RCAT()

    fontsize = 20

    fig = plt.figure(figsize=(8,8))
    ax = plt.axes()
    rcat = R.rcat
    kcat = R.kcat['kcat [s^-1]']
    rmaxn = R.rmaxn['rmax [s^-1]']
    index = kcat.index & rmaxn.index
    x = kcat[index]
    y = rmaxn[index]
    res = np.abs(np.log10(x) - np.log10(y))
    labels = res[res>=1] # at least 10 fold difference
    labels = {k:v for k,v in R.map_reactions_to_gene_names().iteritems() if k in labels}
    report = plot_kcat_rmaxn_correlation(x, y, 
                                         fig, ax,  
                                         labels=labels, 
                                         fit=True)
    
    rmse = np.sqrt( report.sum_square / len(x) )
    r, pval = stats.pearsonr(np.log10(x), np.log10(y))
    
    labels = {k:v for k,v in R.map_reactions_to_gene_names().iteritems() if k in index and k not in labels}    
    add_labels(x, y, labels, ax, fig)    # specific labels to add
    
    ax.set_ylabel(r'in vivo $r_{\mathrm{max}}\,\left[s^{-1}\right]$', 
                  size=fontsize, style='italic')
    ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}\,\left[s^{-1}\right]$', 
                  size=fontsize, style='italic')
    ax.tick_params(axis='both', which='both', top='off', right='off')
    
    [tick.label.set_fontsize(fontsize) for tick in ax.xaxis.get_major_ticks()]
    [tick.label.set_fontsize(fontsize) for tick in ax.yaxis.get_major_ticks()]

    ax.set_xlim(1e-3/4,4*1e3)
    ax.set_ylim(1e-3/4,4*1e3)
    
    plt.tight_layout()
    plt.savefig('%s/svg/kcat_rmax_correlation.svg'%R.path)