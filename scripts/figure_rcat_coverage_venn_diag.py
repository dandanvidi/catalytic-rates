from rcat import RCAT
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3
import csv
from collections import Counter
import numpy as np
from collections import defaultdict

R = RCAT()

genes = R.p
genes = genes[genes>=10]
expressed_genes = set(genes.dropna(how='all').index)

enzymes = set(map(lambda x: x.id, R.model.genes))
enzymes.discard('s0001')
expressed_enzymes = enzymes.intersection(expressed_genes)

homomeric_enzymes = set(R.reactions_by_homomeric_enzymes().values())
homomeric_enzymes.discard('s0001')

homomeric_that_support_flux = defaultdict(list)
for r, b in R.reactions_by_homomeric_enzymes().iteritems():
    if r in R.v.index:
        if R.v.loc[r].sum() > 0:
            homomeric_that_support_flux[b].append(r)
homomeric_that_support_flux = set(homomeric_that_support_flux.keys())
homomeric_that_support_flux.discard('s0001')

fig = plt.figure()
ax = plt.axes()
venn3([enzymes, expressed_enzymes, homomeric_enzymes], 
      ['genes of enzymes', 'expressed enzymatic genes', 'homomeric enzymes'], 
        set_colors=('g', 'b', '#9f4d08'),
        ax=ax)
plt.tight_layout()

not_explained = homomeric_enzymes & expressed_enzymes-homomeric_that_support_flux

# first hypothesis - not explained enzymes support flux in reactions
# in which they do not operate as homomers, i.e., as part of a complex 
# or they catalyze several metabolic reactions

 
#genes = map(model.genes.get_by_id, no_flux_as_homo)
#reactions = map(lambda x: x.reactions, genes)
#no_flux_as_homo_dict = dict(zip(no_flux_as_homo, reactions))
#
#reactions_flux = pd.DataFrame.from_csv('../cache/pFBA_for_all_reactions.csv')
#reactions_flux = reactions_flux[reactions_flux>0].dropna(how='all')
#
## check if some of the expressed homomeric enzymes support flux in other 
## metabolic reactions where they are not unique
#support_flux_as_not_homo = set([])
#for b,rxns in no_flux_as_homo_dict.iteritems():
#    rids = set(map(lambda r: r.id, rxns))
#    test = rids & set(reactions_flux.index)
#    if test != set([]):
#        support_flux_as_not_homo.add(b)
#

#def bnumber_to_functional_group(bnumbers_list=[]):
#
#    b_to_ko_mapping = pd.DataFrame.from_csv('../data/eco_mapping.csv', sep='\t').ko
#    ko_gene_hierarchy_fname = '../data/KO_gene_hierarchy.tms'    
#    systematic_level = 4
#    k_to_function = {}    
#    for row in csv.reader(open(ko_gene_hierarchy_fname, 'r'), delimiter='\t'):
#        if len(row) < systematic_level:
#            continue
#        elif len(row) == systematic_level:
#            function = row[-1]
#        elif len(row) > systematic_level:
#            k_to_function[row[-1]] = function
#    
#    b_to_function = {}    
#    for b in bnumbers_list:
#        try:    
#            b_to_function[b] = k_to_function[list(set((b_to_ko_mapping[b].values)))[0]]
#        except:
#            b_to_function[b] = 'None'
#        
#    return b_to_function
#    
#not_explained = no_flux_as_homo-support_flux_as_not_homo
#b_to_function = bnumber_to_functional_group(not_explained)    
#new_explained = set([x for x in not_explained 
#                            if b_to_function[x] in 
#                            ['tRNA loading', 
#                            'Insulin signaling pathway', 
#                            'Pertussis', 
#                            'Chaperones and folding catalysts',
#                            'Two-component system',
#                            'Phosphatidylinositol signaling system',
#                            'Sulfur relay system']])    
#    
#expressed_enzymes -= support_flux_as_not_homo.union(new_explained)
#
#fig = plt.figure()
#ax = plt.axes()
#venn3([enzymes, expressed_enzymes, homomeric_enzymes], 
#      ['genes of enzymes', 'expressed enzymatic genes', 'homomeric enzymes'], 
#        set_colors=('g', 'b', '#9f4d08'),
#        ax=ax)
#plt.tight_layout()
#plt.savefig('../res/venn_diagram_for_rcat1.svg')
#
#
#fig = plt.figure()
#ax = plt.axes()
#venn3([support_flux, expressed_enzymes, homomeric_enzymes], 
#      ['support flux', 'expressed enzymatic genes', 'homomeric enzymes'], 
#        set_colors=('g', 'b', '#9f4d08'),
#        ax=ax)
#plt.tight_layout()
#plt.savefig('../res/venn_diagram_for_rcat2.svg')
#
#
#
#
#
#import sys, os
#sys.path.append(os.path.expanduser('~/git/across_projects'))
#from plot_types import cdf
#
#kcat = R.get_kcat_of_model_reactions()
#rmax = R.get_rcat_max(7)
#reactions = rmax.index & kcat.index
#
#x = kcat[reactions]
#y = rmax[reactions]
#
#convert_to_irreversible(model)
#residual = np.log10(x/y)
#reac = map(model.reactions.get_by_id, residual.index)
#gens = map(lambda x: list(x.genes)[0], reac)
#revert_to_reversible(model)
#number_of_reactions = map(lambda g: len(g.reactions), gens)
#more_than_one = [r.id for i, r in enumerate(reac) if number_of_reactions[i]>1]
#one_reaction = [r.id for i, r in enumerate(reac) if number_of_reactions[i]==1]
#
#plt.figure()
#ax = plt.axes()
#
#cdf(residual[more_than_one], color='b', ax=ax)
#cdf(residual[one_reaction], color='r', ax=ax)

#
#
#
## highly expressed not explained
#B = E.set_index('bnumber')
#B = B.loc[no_flux_as_homo-support_flux_as_not_homo]
#B.drop('gene_name', axis=1, inplace=True)
#B.drop_duplicates(inplace=True)
#B = B.max(axis=1)
#B.to_csv('../res/expressed_enzymes_with_no_flux.tsv', sep='\t')
##


#
#
