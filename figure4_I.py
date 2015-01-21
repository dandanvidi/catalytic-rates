from kapp import *
import matplotlib.pyplot as plt

log_ticklabels = ['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$10^{0}$',
                    '$10^1$', '$10^2$', '$10^3$', '$10^4$']

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
rate = RCAT(model)
p = PLOT_DATA(model)

kcat = rate.get_kcat_of_model_reactions()

sorted_rates = rate.get_sorted_rates()


reactions = kcat.index & rate.calculate_enzyme_rates().index 
reactions = rate.manually_remove_reactions(reactions)

kcat = kcat[reactions]

pfva_ranges = pd.DataFrame.from_csv('cache/pFVA_ranges_at_1_000000.csv', 
            index_col=[0]) * 6.02214129e23 / 1000 / 3600 #convert units from mmol/gCDW/h to molecules/gCDW/s

idx = pd.MultiIndex.from_tuples([(x.split('\'')[1], x.split('\'')[3]) for x in pfva_ranges.index])

pfva_ranges.index = idx

kcat = kcat[reactions]
rcat_max = rcat_max[reactions]

rmax = pd.Series()
std = pd.Series()
for r in reactions:
    E = rate.E_data.loc[r]
    rcat = pfva_ranges.loc[r].div(E)
        
    rmax[r] = rcat.mean().max()
    
    std[r] = np.std(np.log10(np.sort(rcat.mean().dropna().values)[-3:]))

fig = plt.figure(figsize=(6,6))
ax = plt.axes(axisbg='0.95')
p.plot_kcat_rcat_correlation(x=kcat, y=rmax, fig=fig, ax=ax, color='#FFB84D',
                           edge='none', yerr=10**std, 
                                              labels=reactions)

ax.set_ylabel(r'in vivo $r_{\mathrm{cat, \/max}}$ $\left[s^{-1}\right]$', size=20, style='italic')
ax.set_xlabel(r'in vitro $k_{\mathrm{cat}}$ $\left[s^{-1}\right]$', size=20)
ax.tick_params(axis='both', which='both', top='off', right='off')
ax.set_xticklabels(log_ticklabels, size=12)
ax.set_yticklabels(log_ticklabels, size=12)
#
plt.tight_layout()
plt.show()

fig.savefig('res/figures/kcat_to_rmax_true_max_std_log.pdf')
