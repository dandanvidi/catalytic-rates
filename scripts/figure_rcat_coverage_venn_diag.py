from kapp import RCAT
from cobra.io.sbml import create_cobra_model_from_sbml_file
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3
from cobra.manipulation.modify import convert_to_irreversible
import csv
from collections import Counter
import numpy as np
from cobra.manipulation.modify import revert_to_reversible

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)


hein = pd.DataFrame.from_csv("../data/Heinmann_proteomics.csv") # copies/cell
vilu = pd.DataFrame.from_csv("../data/Vilu_proteomics.csv") # copies/fl[cell]

growth_conditions = pd.DataFrame.from_csv("../data/growth_conditions.csv")

cell_volumes = growth_conditions.cell_volume_l # cell volume in liter
correction = 2 # the data in the table of cell volumes is biased by a factor of 2
hein = hein / ((cell_volumes[hein.columns]/correction) * 1e15) # convert from # copies/cell to copies/fl[cell]

expression = hein.join(vilu, how='outer')
expression = expression[growth_conditions.index]

''' now both data sets are in copies/fl cytoplasm '''
expression = (expression[expression>10])# remove genes with less then 10 copies per fl cytoplasm
expression.dropna(how='all', inplace=True)

expressed_genes = set(expression.index)

R = RCAT(model)
E = R.E_data
V = R.V_data

enzymes = set(map(lambda x: x.id, model.genes))
homomeric_enzymes =set(E.bnumber) #668
expressed_enzymes = enzymes & expressed_genes
homomeric_enzymes.discard('s0001')
enzymes.discard('s0001')

support_flux = set([])
for b in homomeric_enzymes:
    reactions = V[V.bnumber==b]
    if reactions.sum(axis=1).sum() > 0: #at least one reaction supports flux:
        support_flux.add(b)
#support_flux.discard('s0001')

no_flux_as_homo = homomeric_enzymes & expressed_enzymes-support_flux
genes = map(model.genes.get_by_id, no_flux_as_homo)
reactions = map(lambda x: x.reactions, genes)
no_flux_as_homo_dict = dict(zip(no_flux_as_homo, reactions))

reactions_flux = pd.DataFrame.from_csv('../cache/pFBA_for_all_reactions.csv')
reactions_flux = reactions_flux[reactions_flux>0].dropna(how='all')

# check if some of the expressed homomeric enzymes support flux in other 
# metabolic reactions where they are not unique
support_flux_as_not_homo = set([])
for b,rxns in no_flux_as_homo_dict.iteritems():
    rids = set(map(lambda r: r.id, rxns))
    test = rids & set(reactions_flux.index)
    if test != set([]):
        support_flux_as_not_homo.add(b)


def bnumber_to_functional_group(bnumbers_list=[]):

    b_to_ko_mapping = pd.DataFrame.from_csv('../data/eco_mapping.csv', sep='\t').ko
    ko_gene_hierarchy_fname = '../data/KO_gene_hierarchy.tms'    
    systematic_level = 4
    k_to_function = {}    
    for row in csv.reader(open(ko_gene_hierarchy_fname, 'r'), delimiter='\t'):
        if len(row) < systematic_level:
            continue
        elif len(row) == systematic_level:
            function = row[-1]
        elif len(row) > systematic_level:
            k_to_function[row[-1]] = function
    
    b_to_function = {}    
    for b in bnumbers_list:
        try:    
            b_to_function[b] = k_to_function[list(set((b_to_ko_mapping[b].values)))[0]]
        except:
            b_to_function[b] = 'None'
        
    return b_to_function
    
not_explained = no_flux_as_homo-support_flux_as_not_homo
b_to_function = bnumber_to_functional_group(not_explained)    
new_explained = set([x for x in not_explained 
                            if b_to_function[x] in 
                            ['tRNA loading', 
                            'Insulin signaling pathway', 
                            'Pertussis', 
                            'Chaperones and folding catalysts',
                            'Two-component system',
                            'Phosphatidylinositol signaling system',
                            'Sulfur relay system']])    
    
expressed_enzymes -= support_flux_as_not_homo.union(new_explained)

fig = plt.figure()
ax = plt.axes()
venn3([enzymes, expressed_enzymes, homomeric_enzymes], 
      ['genes of enzymes', 'expressed enzymatic genes', 'homomeric enzymes'], 
        set_colors=('g', 'b', '#9f4d08'),
        ax=ax)
plt.tight_layout()
plt.savefig('../res/venn_diagram_for_rcat1.svg')


fig = plt.figure()
ax = plt.axes()
venn3([support_flux, expressed_enzymes, homomeric_enzymes], 
      ['support flux', 'expressed enzymatic genes', 'homomeric enzymes'], 
        set_colors=('g', 'b', '#9f4d08'),
        ax=ax)
plt.tight_layout()
plt.savefig('../res/venn_diagram_for_rcat2.svg')





import sys, os
sys.path.append(os.path.expanduser('~/git/across_projects'))
from plot_types import cdf

kcat = R.get_kcat_of_model_reactions()
rmax = R.get_rcat_max(7)
reactions = rmax.index & kcat.index

x = kcat[reactions]
y = rmax[reactions]

convert_to_irreversible(model)
residual = np.log10(x/y)
reac = map(model.reactions.get_by_id, residual.index)
gens = map(lambda x: list(x.genes)[0], reac)
revert_to_reversible(model)
number_of_reactions = map(lambda g: len(g.reactions), gens)
more_than_one = [r.id for i, r in enumerate(reac) if number_of_reactions[i]>1]
one_reaction = [r.id for i, r in enumerate(reac) if number_of_reactions[i]==1]

plt.figure()
ax = plt.axes()

cdf(residual[more_than_one], color='b', ax=ax)
cdf(residual[one_reaction], color='r', ax=ax)

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
