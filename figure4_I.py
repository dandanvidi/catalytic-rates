''' plots the correlation plot with standard deviations as the mean rcat and max rcat'''

from kapp import *
import matplotlib.pyplot as plt

log_ticklabels = ['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$',
                    '$10^1$', '$10^2$', '$10^3$', '$10^4$']

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
rate = RCAT(model)
p = PLOT_DATA(model)

kcat = rate.get_kcat_of_model_reactions()
rcat = rate.calculate_enzyme_rates()    
rcat_max = rate.get_rcat_second_max()
rcat_median = rcat.median(axis=1)    

sorted_rates = rate.get_sorted_rates()


reactions = kcat.index & rcat_max.index 
reactions = rate.manually_remove_reactions(reactions)

pfva_ranges = pd.DataFrame.from_csv('cahce/pFVA_ranges_at_1.010.csv')

flux_uncertainty = {}
for r in reactions:
    cond = rate.get_best_condition()[r]
    flux_uncertainty[r] = pfva_ranges[cond][]

res = rcat_max[reactions] / kcat[reactions]

yerr = pd.Series({k:flux_uncertainty[c][k] * rcat_max[k]
                    for k,c in rate.get_best_condition().iteritems()
                    if k in reactions})
                        
   
#below = np.log10(pd.Series({k:flux_uncertainty[c][k]*100*rcat_max[k]/2  
#                    for k,c in rate.get_best_condition().iteritems()
#                    if k in reactions}))
                        

fig = plt.figure(figsize=(6,6))
ax = plt.axes(axisbg='0.95')
p.plot_kcat_rcat_correlation(x=kcat[r], y=rcat_max[r], fig=fig, ax=ax, color='#FFB84D',
                           edge='none', yerr=yerr, labels=reactions)

ax.set_ylabel(r'in vivo $r_{\mathrm{cat, \/max}}$ $\left[s^{-1}\right]$', size=20, style='italic')
ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}$ $\left[s^{-1}\right]$', size=20)
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.set_xticklabels(log_ticklabels, size=12)
ax.set_yticklabels(log_ticklabels, size=12)

plt.tight_layout()
plt.show()

fig.savefig('res/kcat_to_rcat_max_correlation.pdf')
