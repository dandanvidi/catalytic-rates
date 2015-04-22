from cobra.core import Metabolite, Reaction
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.variability import flux_variability_analysis
import os
from scipy import stats, odr
import math
import json
import csv
import sys
from collections import defaultdict
from xml.dom.minidom import parse
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import numpy as np
import pandas as pd
from scipy.stats.mstats import gmean
sys.path.append(os.path.expanduser('~/git/component-contribution'))
from python.kegg_reaction import KeggReaction
from python.kegg_model import KeggModel
from python.component_contribution import ComponentContribution
from python.thermodynamic_constants import R, default_T
from copy import deepcopy

class MODEL(object):
    
    def __init__(self, model):
        """
        model - a cobra object
        model_fname - str of path to model xml file
        """
        self.model = deepcopy(model)
        self.gene_names = {row[0:5]:row[84:].split(';')[0].strip() 
                            for row in open("../data/all_ecoli_genes.txt", 'r')}
        
    def model_metabolites(self):
        
        '''map model metabolites to kegg CIDs'''
        metabolites = pd.read_csv("../data/metaboliteList.txt", delimiter='\t')
        metabolites.pop('cas_id')
        metabolites.pop('Unnamed: 8')
        metabolites.set_index('name', inplace=True)
        return metabolites

    def reaction_formula(self, reaction_list):
        '''
            reaction_list - list of cobra reaction objects
            returns a dictionary
                keys: reactions
                value: reactions str in kegg CIDs
        '''

        reaction_list = map(self.model.reactions.get_by_id, reaction_list)
        CIDS = self.model_metabolites().kegg_id.dropna().astype('str')

#        #replace all reactions sparse with CID sparses
        r_dict = {}
        sparse = {}
        reaction_strings = {}
        for r in reaction_list:
            r_dict = {CIDS[m.name]:float(r.metabolites[m]) 
                            for m in r.metabolites
                            if m.name in CIDS.index}
            if len(r_dict) != len(r.metabolites):
                continue

            if 'C00080' in r_dict:
                del r_dict['C00080']            
            sparse[r.id] = r_dict
            assert r not in reaction_strings, 'Duplicate reaction!'
            kegg_reaction = KeggReaction(r_dict)
#            print kegg_reaction

#            kegg_reaction.remove_protons()
            reaction_strings[r.id] = str(kegg_reaction)

        return reaction_strings, sparse

    def map_model_reactions_to_EC(self, model_fname):
        
        ''' parse the sbml model file to extract EC numbers of reactions '''
        
        document = parse(model_fname)
        ec_list = []
        for r_elem in document.getElementsByTagName('reaction'):
            for p_elem in r_elem.getElementsByTagName('p'):
                val = p_elem.childNodes[0].nodeValue
                if val.find('EC Number:') != -1:
                    ec = val[11:].encode('utf-8')            
                    ec_list.append(ec)

        r_to_ec = dict(zip(self.model.reactions, ec_list))        
        convert_to_irreversible(self.model)
        for r in self.model.reactions:
            if 'reverse' in r.id:
                original_r = self.model.reactions.get_by_id(r.id[:-8])
                r_to_ec[r] = r_to_ec[original_r]
                
        r_to_ec = {r.id:v for r,v in r_to_ec.iteritems()}    
        r_to_ec = pd.DataFrame(r_to_ec.items()).set_index(0)
        r_to_ec.to_csv("../cache/reactions_to_ec_mapping.csv")
        
        return 
       
    def map_model_reaction_to_genes(self):
        ''' 
            find and match all reactions in the model 
            that are catalyzed by a single gene 
        '''

        convert_to_irreversible(self.model)
        r_to_single_gene = pd.Series()

        # remove reactions which can be catalysed by more than one gene
        enzymatic_reactions = 0
        for r in self.model.reactions:
            genes = map(lambda x: x.id, r.genes)
            if len(genes) > 0:
                enzymatic_reactions += 1
            if len(genes) == 1 and genes[0]:
                r_to_single_gene[r.id] = genes[0]
        
        # manual additions of genes with only one active isoenzyme 
        r_to_single_gene["METS"] = "b3829" # metE - cobalamin-independent homocysteine transmethylase
        r_to_single_gene["HCO3E"] = "b0126" # can - carbonic anhydrase
        r_to_single_gene["PFK"] = "b3916" #  6-phosphofructokinase                
        r_to_single_gene["RPI"] = "b2914" #  ribose-5-phosphate isomerase A                
        
#        r_to_single_gene.to_csv("cache/reactions_to_genes_mapping.csv")
#
        return r_to_single_gene
#        
    def map_genes_to_model_reaction(self):
        ''' 
            find and match all genes in E. coli to reactions in the model 
            that are catalyzed by a single gene
        '''

        convert_to_irreversible(self.model)
        r_to_single_gene = pd.Series()

        # remove reactions which can be catalysed by more than one gene
        for r in self.model.reactions:
            genes = map(lambda x: x.id, r.genes)
            if len(genes) == 1 and genes[0]:
                r_to_single_gene[genes[0]] = r.id
        
        # manual additions of genes with only one active isoenzyme 
        r_to_single_gene["b3829"] = "METS" # metE - cobalamin-independent homocysteine transmethylase
        r_to_single_gene["b0126"] = "HCO3E"# can - carbonic anhydrase
        r_to_single_gene["b3916"] = "PFK" #  6-phosphofructokinase                
        r_to_single_gene["b2914"] = "RPI" #  ribose-5-phosphate isomerase A                
        
#        r_to_single_gene.to_csv("cache/genes_to_reactions_mapping.csv")
#
        return r_to_single_gene
        
class RCAT(MODEL):
    
    def __init__(self, model):
        MODEL.__init__(self, model)

        self.growth_conditions = pd.DataFrame.from_csv("../data/growth_conditions.csv")
        self.growth_conditions.sort(axis=0, columns='growth_rate_1_h', inplace=True)

        self.growth_conditions.drop('anaerobic_heinmann', axis=0, inplace=True)
        
        self.V_data = pd.DataFrame.from_csv('../cache/fluxes_[molecules_per_second_per_gCDW].csv')        
        self.E_data = pd.DataFrame.from_csv("../cache/expression_[copies_per_gCDW].csv")
        
        columns = ['bnumber', 'gene_name'] + list(self.growth_conditions.index)
        
        self.V_data = self.V_data[columns]
        self.E_data = self.E_data[columns]
        
    def _get_enzyme_abundance(self, proteomics_fname):
        '''
            map proteomics data to the cobra model reactions
            by bnumbers (gene identifiers in E. coli)
            
            Argumetns:
                csv file with proteomics data
                units are: copies of enzmye / gCDW
            
            Returns:
                pandas series with reactions as index and
                enzyme abundance as values
        '''

        reactions_to_bnumbers = self.map_model_reaction_to_genes().set_index(1)
        new_data = pd.read_csv(proteomics_fname).set_index('bnumber')
        
        abundance = reactions_to_bnumbers.join(new_data, how='left')
        abundance.set_index(0, inplace=True)
        abundance.dropna(inplace=True)
        
        return abundance

    def _calculate_pFBA(self, title, growth_params):
        '''
            calculates parsimoniuos FBA - the objective function is
            to minimize the total sum of fluxes
            
            Arguments:
                growth parametres:
                    carbon source
                    growth_rate [1/h]
                    oxygen uptake rate [mmol/gCDW/h]
                    
            Returns:
                pandas DataFrame flux distribution of COBRA model reactions
                in units of molecules/gCDW/s
                
        '''        

        convert_to_irreversible(self.model)            

        cond = growth_params['carbon']
        growth_rate = growth_params['growth_rate_1_h']
        oxy = growth_params['oxygen_uptake_mmol_gCDW_h']
                
        rxns = dict([(r.id, r) for r in self.model.reactions])
        rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
        try:
            rxns['EX_' + cond + '_e'].lower_bound = -18.5 # redefine sole carbon source uptake reaction in mmol/gr/h
        except:
            rxns['EX_glc_e'].lower_bound = -1000
        rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = growth_rate            
        rxns['EX_o2_e'].lower_bound = oxy
        print "solving model for %s..." %title,
        optimize_minimal_flux(self.model, already_irreversible=True)
        print "\t GR: %.02f 1/h" %self.model.solution.f
        
        flux_dist = pd.DataFrame(self.model.solution.x_dict.items()).set_index(0)
        
        #convert units from mmol/gCDW/h to molecules/gCDW/s
        flux_dist = flux_dist * 6.02214129e23 / 1000 / 3600 

        return flux_dist    
        
    def calculate_enzyme_rates(self):
        '''
            calculates the catalytic rate of enzymes by dividing the flux by
            the copy number of enzymes. 
            
            Important: 
                Data must have matching units
        '''
        V = self.V_data.drop(['bnumber', 'gene_name'], axis=1)        
        E = self.E_data.drop(['bnumber', 'gene_name'], axis=1)
        
        # calculate rats by dividing flux by proteomics
        rcat = V / E

        # replace zero and inf values with nan
        rcat.replace([0, np.inf, -np.inf], np.nan, inplace=True)
        rcat.dropna(how='all', inplace=True)

        return rcat

    def get_rcat_max(self, minimal_conditions=5):

        rcat = self.calculate_enzyme_rates()
        rcat.dropna(thresh=minimal_conditions, inplace=True)
        return rcat.max(axis=1)
        
    def get_cond_max(self, minimal_conditions=5):

        rcat = self.calculate_enzyme_rates()
        rcat.dropna(thresh=minimal_conditions, inplace=True)
        return rcat.idxmax(axis=1)

    def get_sorted_rates(self):
        '''
            sorts the catalytic rates of all (enzyme, reaction) pairs
            
            Returns:
                Dictionary woth reaction as keys and sorted array as value
        '''
        rcat = self.calculate_enzyme_rates()
        
        rcat = rcat.T
        sorted_rates = {}
        for r in rcat:
            n = rcat[r].dropna().values
            sorted_rates[r] = sorted(n, reverse=True)

        return sorted_rates
        
    def get_kcat_of_model_reactions(self):
        
        # kcat values collected from BRENDA and other publications - manually curated
        kcats = pd.read_csv("../data/manual_kcat_values.csv")
        kcats.set_index('reactions', inplace=True)

        return kcats['kcat per subunit']        
        
    def add_condition(self, proteomics_fname, growth_params):
        '''
            growth params  - a dictionary of growth condition paramters
                title - name of condition
                carbon - str of carbon source in a cobra style (e.g glucose -> glc)
                growth_rate - float of growht rate [1/h]
                oxygen - float of oxygen uptake rate [mmol/gCDW/h]
        '''
        
        new_abundance = self._get_enzyme_abundance(proteomics_fname)
        out_abundance = self.E_data.join(new_abundance, how='outer')
        out_abundance.to_csv("cache/enzyme_conc_across_conditions.csv")
        
        new_flux = self._calculate_pFBA(growth_params)
        out_flux = self.V_data.join(new_flux, how='outer')
        out_flux.to_csv("cache/pFBA_distribution_across_conditions.csv")
                
    def delete_non_expressed_genes(self):
        
        """
            delete genes from model if they are not represented in the 
            proteomics data of a given condition
            
            A
        """
    
    def calculate_pFVA(self, title, growth_params, relaxation=1):
        '''
            calculates the uncertainties of flux by running FVA 
            given minimization of the total sum of fluxes
            
            Arguments:
                growth_params - 
                
            Returns:
                flux_ranges - the range of each flux
        '''
        
        model = deepcopy(self.model)        
        cond = growth_params['carbon']
        growth_rate = growth_params['growth_rate_1_h']
        oxy = growth_params['oxygen_uptake_mmol_gCDW_h']
        
        convert_to_irreversible(model)

        rxns = dict([(r.id, r) for r in model.reactions])

        rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
        rxns['EX_' + cond + '_e'].lower_bound = -1000 # redefine sole carbon source uptake reaction in mmol/gr/h

        rxns['EX_o2_e'].lower_bound = oxy

        rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = growth_rate            
        rxns['Ec_biomass_iJO1366_core_53p95M'].lower_bound = growth_rate            

        fake = Metabolite(id='fake')
        model.add_metabolites(fake)        
                
        for r in model.reactions:
            r.add_metabolites({fake:1})
            
        flux_counter = Reaction(name='flux_counter')
        flux_counter.add_metabolites(metabolites={fake:-1})                

        model.add_reaction(flux_counter) 
        model.change_objective(flux_counter)
        
        print "solving pFVA -> " + title + " -> %f relaxation" %relaxation
        fva = flux_variability_analysis(model, 
                                        reaction_list=model.reactions, 
                                        objective_sense='minimize',
                                        fraction_of_optimum=relaxation)
        return fva

    def manually_remove_reactions(self, reactions):
        '''
            some reactions are mannualy removed from the analyis
            as there might be issues with annotations 
            or error prone flux calculation
            
            Arguments:
                reaction list
                
            Returns:
                new list of reactions with the problematic reactions removed
        '''
                
        # PPKr_reverse reaction is used for ATP generation from ADP 
        # in the FBA model. Nevertheless, acording to EcoCyc, it is used to 
        # to generate polyP (inorganic phosphate) chains from ATP and it is not
        # part of the oxidative phosphorilation
        list_of_reactions = ['PPKr_reverse']
        return reactions.drop(list_of_reactions) 
        
    def calculate_rcat_err(self):
        s = self.get_sorted_rates()
        s = {k:np.std(v[0:2]) for k, v in s.iteritems()}
        s = pd.DataFrame(s, index=['std'])
        s.index.name = 'reactions'
        return s.T

class MM_KINETICS(RCAT):

    def __init__(self, model):
        RCAT.__init__(self, model)    
        self.CC_CACHE_FNAME = os.path.expanduser('../../component-contribution/cache/component_contribution.mat')

        self.reactions = [r.id for r in self.model.reactions]
        self.metab_conc = pd.read_csv('../data/metab_conc.csv') 
        self.metab_conc.set_index('Compound Id (KEGG)', inplace=True)
        self.metab_conc.dropna(how='all', inplace=True)
        
        self.known_cids = {'C'+'%05i'%c:c for c in self.metab_conc.index}
        self.KMs = pd.DataFrame.from_csv('../cache/KM_values.csv')

    def get_reactions_metabolites(self):
        self.r_to_metab = pd.DataFrame(columns=['substrates', 'products'],
               index=self.reactions)
                                       

        for r in self.model.reactions:
            substrates = [m.name for m, i in r.metabolites.iteritems() 
                                    if i < 0]                                       
            products = [m.name for m, i in r.metabolites.iteritems() 
                                    if i > 0]
            self.r_to_metab['substrates'][r.id] = substrates
            self.r_to_metab['products'][r.id] = products

    def cache_reactions_dG(self, c = 1e-3, pH=7.5, I=0.2, T=298.15):
    
        if not os.path.exists(self.CC_CACHE_FNAME):
            cc = ComponentContribution()
            cc.save_matfile(self.CC_CACHE_FNAME)
        else:
            cc = ComponentContribution.from_matfile(self.CC_CACHE_FNAME)    
        
        convert_to_irreversible(self.model)   
        model_reactions = [r.id for r in self.model.reactions]
        reaction_strings, sparse = self.reaction_formula(model_reactions)

        reactions = []
        reac_strings = []
        for key, val in reaction_strings.iteritems():    
            reactions.append(key)
            reac_strings.append(val)
            
        Kmodel = KeggModel.from_formulas(reac_strings)
        Kmodel.add_thermo(cc)
        dG0_prime, dG0_std = Kmodel.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
        # use 1mM as the standard concentration
        print c
        conc = np.ones((1, len(Kmodel.cids))) * 1e-3

        for cid in ['C00001']: # set the conc of water as 1M
            i = Kmodel.cids.index(cid)
            conc[0, i] = 1
            
        r_to_dGc = pd.DataFrame(index = reactions, columns = ['glc', 'glyc', 'ac', 'std'])
        for condition in ['glc', 'glyc', 'ac']:
            for ID, cid  in self.known_cids.iteritems(): # set the conc of known metabolites
                if ID in Kmodel.cids:
                    i = Kmodel.cids.index(ID)
                    c = self.metab_conc[condition][cid]
                    if not np.isnan(c):
                        conc[0, i] = c
                
            dGc_prime = dG0_prime + R * default_T * np.dot(np.log(conc), Kmodel.S).T
            r_to_dGc[condition] = dGc_prime

#        print len(dG0_std)
#        r_to_dGc['std'] = dG0_std
        r_to_dGc.to_csv('../cache/reactions_to_dGc.csv')

    def _backwards_reaction_effect(self):

        dG = pd.DataFrame.from_csv('../cache/reactions_to_dGc.csv')
        dG = dG[['glc','glyc','ac']]
        T = 1-np.exp(dG/(R*default_T))
        T.replace([np.inf, -np.inf], np.nan, inplace=True)

        return dG, T
        
    def _calculate_saturation_effect(self, kms, stoicho, condition):
        
        metabolites = kms.index
        smul = 1
        for m in metabolites:
            cid = int(self.model_metabolites()['kegg_id'][m][1:])
    #        conc = 100 # uM
            if cid in self.metab_conc[condition].dropna():
                conc = self.metab_conc[condition][cid]*1000 
                s = (conc / kms[m])**stoicho[m]
                smul *= s
                return (smul) / (1 + smul)


    def concentraion_dependant_effects(self):
        kcat = self.get_kcat_of_model_reactions()
        rmax = self.get_rcat_max(7)   
        relevant_reac = set(rmax.index & kcat.index & self.KMs.index)
        cmax = self.get_cond_max(7)
        convert_to_irreversible(self.model)
        carbon_cond = {r:self.growth_conditions['carbon'][cmax[r]] for r in cmax.index}
        dG, thermo = self._backwards_reaction_effect()
        
        S = pd.Series()
        T = pd.Series()
        for r in self.model.reactions:
            if r.id in relevant_reac:
                stoicho = {k.name:-v for k,v in r.metabolites.iteritems() 
                                                if v<0 and k.name not in ['H2O', 'H+']}
                kms = self.KMs.loc[r.id].dropna()
                kms = kms[kms<0]*-1
                condition = carbon_cond[r.id]
                if condition in ['glc', 'ac', 'glyc']:
                    if len(kms) == len(stoicho):
                        S[r.id] = self._calculate_saturation_effect(kms, stoicho, condition)
                    if r.id in thermo.index:
                        T[r.id] = thermo.loc[r.id, condition]        
        
        T = np.abs(T[T>-0.1])
        S.dropna(inplace=True)
        S = S.astype('float')
        
        intersect = T.index & S.index
        S = S[intersect]
        T = T[T>0][intersect]
        ST = (S*T).dropna()
        
        return S, T, ST
        
class PLOT_DATA(RCAT):

    def __init__(self, model):
        RCAT.__init__(self, model) 

    def expression_plot(self, ax, e_change, v_change=0, legend=False):
    
        ax.hist(e_change.dropna(), bins=np.arange(-3, 4.01, 0.5), zorder=9, label=r'$\frac{E_n}{E_m}$', 
                facecolor=(1.0,0.6,0.6), edgecolor=(0.5, 0.3, 0.3), histtype='stepfilled')
        ax.axvline(e_change.median(), 0, 1, c=(0.6, 0, 0), zorder=10, ls=':', 
                   lw=3, label=r'$\frac{E_n}{E_m}$ $\rm{median}$')
#        ax.axvline(v_change.median(), 0, 1, zorder=10, lw=3, c='k', ls='-', 
#                   label=r'$\frac{v_n}{v_m}$ $\rm{median}$')
        
        if legend:
            ax.legend(fontsize=15, bbox_to_anchor=(0.175, 0.8))
        
        xlabels = [2**i for i in np.arange(-3,0)] + [2**i for i in np.arange(0,5)]
        ax.set_xlim(-3,4)
        ax.set_xticklabels(xlabels)
#        ax.arrow(e_change.median(), 30, v_change.median()-e_change.median(), 0,
#                 head_width=0.05, head_length=0.1,
#                 fc='k', ec='k')


    def plot_kcat_rcat_correlation(self, x, y, fig, ax, color='b', edge='none', 
                                   yerr='none', labels=[], fit_on=True):
    
        logx = np.log10(x)
        logy = np.log10(y)
        
        ax.scatter(x, y,s=40, c=color, alpha=0.7, marker='o', edgecolor=edge, zorder=3)
        
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
    
        if fit_on:  
            edge = np.array([-4, 4])
            ax.plot([1e-4, 1e4], 10**fit_func(beta, edge), color='#699A33', ls=':', lw=3, zorder=1)
            
            
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        b_to_r = self.map_model_reaction_to_genes()
        
        if labels!=[]:
            ann = []
            for r in labels:
                if x[r]>y[r]:
                    ann.append(ax.text(x[r]*1.2, y[r], 
                                       self.gene_names[b_to_r[r]], 
                                        ha='left', va='center', zorder=5))
                if x[r]<y[r]:
                    ann.append(ax.text(x[r]/1.3, y[r], 
                                       self.gene_names[b_to_r[r]], 
                                        ha='right', va='center', zorder=5))
                                        
            mask = np.zeros(fig.canvas.get_width_height(), bool)
            fig.canvas.draw()
            for i, a in enumerate(ann):
                bbox = a.get_window_extent()
                x0 = int(bbox.x0)
                x1 = int(math.ceil(bbox.x1))
                y0 = int(bbox.y0)
                y1 = int(math.ceil(bbox.y1))
            
                s = np.s_[x0:x1, y0:y1]
                if np.any(mask[s]):
                    a.set_visible(False)
                else:
                    mask[s] = True
                    
        for tickx, ticky in zip(ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()):
            tickx.label.set_fontsize(14) 
            ticky.label.set_fontsize(14)
            
        ax.set_xlim(1e-4,1e4)
        ax.set_ylim(1e-4,1e4)
#        ax.set_xticklabels(ax.get_xticklabels(), size=30)
#        ax.grid(color='b', alpha=0.5, ls='-', lw=1, zorder=0)

        cor, pval = stats.pearsonr(logx, logy)
        rsq = cor**2
        ax.text(1e-3/2, 1e3*2, r'$R^2=$%.2f' %rsq, size=15)        
        rmse = np.sqrt( output.sum_square / len(x) )
        print "slope and intercept = " , beta
        print "r^2 = %.3f, pval = %.2e"%(cor**2, pval)
        print "rmse = %.3f" % rmse    
        return output

if __name__ == "__main__":

    model_fname = "../data/iJO1366_curated.xml"
    model = create_cobra_model_from_sbml_file(model_fname)
    mm = MM_KINETICS(model)
    
#    mm.reaction_formula(map(lambda x: x.id, model.reactions))
    o = mm.map_model_reactions_to_EC(model_fname)
#    mm.cache_reactions_dG()

    
