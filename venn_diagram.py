from kapp import RCAT
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import numpy as np
from collections import Counter
import csv

model = create_cobra_model_from_sbml_file("../data/iJO1366_curated.xml")
convert_to_irreversible(model)
R = RCAT(model)

analyzable_reactions = R.V_data

v = R.V_data.drop(['bnumber', 'gene_name'], axis=1)
E = R.E_data.drop(['bnumber', 'gene_name'], axis=1)

support_flux = v[v>0].dropna(how='all').index
not_support_flux = v.index - support_flux

#not_support_flux[not_support_flux != 0] = np.nan
#not_support_flux = not_support_flux.dropna(how='all').index

expressed_enzymes = E.dropna(how='all').index
not_expressed_enzymes = E.index - expressed_enzymes 

FE = support_flux & expressed_enzymes
Fe = support_flux & not_expressed_enzymes
fE = not_support_flux & expressed_enzymes
fe = not_support_flux & not_expressed_enzymes

explained_gene_expression = set(R.V_data['bnumber'].loc[FE].values)
not_explained = []
for i in fE:
    r = model.reactions.get_by_id(i)
    if r.id == 'RPI':
        gene = set(['b2914'])
    else:        
        gene = set(map(lambda x: x.id, r.genes))
    if not gene & explained_gene_expression:
        not_explained.append(r)
        

ne = map(lambda x: x.id, not_explained)

ne_total_fraction = (E.loc[ne].median(axis=1) / E.sum().median()).sum()



genes = []
subsystems = []
for r in not_explained:
    genes += map(lambda x: x.id, r.genes)
    subsystems.append(r.subsystem)
genes = set(genes)

def bnumber_to_functional_group(bnumbers_list=[]):

    b_to_ko_mapping = pd.DataFrame.from_csv('../data/eco_mapping.csv', sep='\t').ko
    ko_gene_hierarchy_fname = '../data/KO_gene_hierarchy.tms'    
    j = 0    
    systematic_level = 3
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
            
b_to_function = bnumber_to_functional_group(genes)    
#eco_mapping = pd.DataFrame.from_csv('../data/eco_mapping.csv', sep='\t')
#i = 0
#for g in genes:
#    try:
#        eco_mapping.ko[g]
#    except KeyError:
#        continue

#b_to_r = R.map_genes_to_model_reaction()
#
#reactions = R.V_data.index
#supported_fluxes = R.V_data.drop(['bnumber', 'gene_name'], axis=1).replace(0, np.nan).dropna(how='all').index
#rcat = R.calculate_enzyme_rates().index
#extra = expressed_enzymes & supported_fluxes - rcat
#
#fig = plt.figure()
#ax = plt.axes()
#venn3([reactions, expressed_enzymes-extra, supported_fluxes], 
#      set_labels=('reactions by single \n homomeric enzyme', 'enzymes expressed', 'flux supported'),
#      ax=ax)
#
#ne = expressed_enzymes -supported_fluxes
#j = 0
#unrecognized_pathways = {}
#unrecognized_genes = {}
#for i in ne:
#    try:
#        r = model.reactions.get_by_id(i)
#        unrecognized_pathways[r] = r.subsystem
#        unrecognized_genes[r] = [R.gene_names[g.id] for g in r.genes]
#    except KeyError:
#        j += 1
#
#genes = []
#for r,glist in unrecognized_genes.iteritems():
#    genes += glist
#gene_count = Counter(genes)

#report = csv.writer(open('../res/unexplained_expression.csv', 'w'), delimiter='\t')
#count = Counter(unrecognized_pathways)
#for r in out.iterkeys():
#    print count[r], unrecognized_pathways
#    report.writerow([s, v])
#    
#plt.tight_layout()
#plt.savefig('../res/venn_diagram.svg')



#hein = pd.DataFrame.from_csv('../data/Heinmann_proteomics.csv')
#vilu = pd.DataFrame.from_csv('../data/Vilu_proteomics.csv')
#genes = set(hein.join(vilu, how='outer').dropna(how='all').index)
#enzymes = set(map(lambda x: x.id, model.genes))
#expressed_enzymes = enzymes & genes
#R = RCAT(model)

#

#
#expressed = set([])
#for b in enzymes & genes:
#    try:
#        expressed.add(b_to_r[b])
#    except KeyError:
#        expressed.add(b)
#        continue
#
#growth_conditions = pd.DataFrame.from_csv('../data/growth_conditions.csv')
#reactions = set([])
#for index in growth_conditions.index:
#    model = create_cobra_model_from_sbml_file("../data/iJO1366_curated.xml")
#    convert_to_irreversible(model)
#    R = RCAT(model)
#    growth_params = growth_conditions.loc[index]
#    flux = R._calculate_pFBA(growth_params)
#    flux = flux[flux>0].dropna()
#    reactions = reactions.union(flux.index)
#    
#
#venn3([enzymes, fluxes.index, rcat.index], set_labels=('expressed', 'reactions', 'rcat'))