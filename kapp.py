from cobra.core import Metabolite, Reaction
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.variability import flux_variability_analysis
import os
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

class MODEL(object):
    
    def __init__(self, model):
        """
        model - a cobra object
        model_fname - str of path to model xml file
        """
        self.model = model
        
    def model_metabolites(self):
        
        '''map model metabolites to kegg CIDs'''
        metabolites = pd.read_csv("data/metaboliteList.txt", delimiter='\t')
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

        reaction_list = map(model.reactions.get_by_id, reaction_list)
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
            
            sparse[r.id] = r_dict
            assert r not in reaction_strings, 'Duplicate reaction!'
            kegg_reaction = KeggReaction(r_dict)
            kegg_reaction.remove_protons()
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
        r_to_ec.to_csv("cache/reactions_to_ec_mapping.csv")
        
        return 
       
    def map_model_reaction_to_genes(self):
        ''' 
            find and match all reactions in the model 
            that are catalyzed by a single gene 
        '''

        convert_to_irreversible(self.model)
        r_to_single_gene = {}

        # remove reactions which can be catalysed by more than one gene
        for r in self.model.reactions:
            genes = map(lambda x: x.id, r.genes)
            if len(genes) == 1 and genes[0]:
                r_to_single_gene[r.id] = genes[0]
        
        # manual additions of genes with only one active isoenzyme 
        r_to_single_gene["METS"] = "b3829" # metE - cobalamin-independent homocysteine transmethylase
        r_to_single_gene["HCO3E"] = "b0126" # can - carbonic anhydrase
        r_to_single_gene["PFK"] = "b3916" #  6-phosphofructokinase                
        r_to_single_gene["RPI"] = "b2914" #  ribose-5-phosphate isomerase A                
        
        r_to_single_gene = pd.DataFrame(r_to_single_gene.items())
        r_to_single_gene.to_csv("cache/reactions_to_genes_mapping.csv")
#
        return r_to_single_gene
#        
class RCAT(MODEL):
    
    def __init__(self, model):
        MODEL.__init__(self, model)

        self.E_data = pd.read_csv("cache/enzyme_conc_across_conditions.csv")
        self.E_data.set_index('reactions', inplace=True)

        self.V_data = pd.read_csv("cache/pFBA_dist_across_conditions.csv")
        self.V_data.set_index('reactions', inplace=True)

    def calculate_enzyme_rates(self):
        
        # find intercept of measured conditions
        if len(self.E_data.columns) != len(self.V_data.columns):
            print "the number of conditions is not identical \
                  in flux and proteomics, thus the intersect is used"

        conditions = self.E_data.columns & self.V_data.columns
        
        # calculate rats by dividing flux by proteomics
        rates = self.V_data[conditions] / self.E_data[conditions]

        # replace zero and inf values with nan
        rates.replace([0, np.inf, -np.inf], np.nan, inplace=True)

        return rates

    def get_rcat_max(self, minimal_conditions):
        
        rcat = self.calculate_enzyme_rates()      
        rcat.dropna(thresh=minimal_conditions, inplace=True)
        return rcat.max(axis=1)       
            
    def get_rcat_second_max(self, minimal_conditions):
        
        rcat = self.calculate_enzyme_rates()
        rcat.dropna(thresh=minimal_conditions, inplace=True)

        max1 = zip(rcat.index, rcat.idxmax(axis=1))
        
        for r, cond in max1:
            rcat[cond][r] = 0
        return rcat.max(axis=1)        

        
    def get_best_condition(self, minimal_conditions):
        
        rcat = self.calculate_enzyme_rates        
        rcat.dropna(thresh=minimal_conditions, inplace=True)
        return rcat.idxmax(axis=1)       
        
    def get_kcat_of_model_reactions(self):
        
        # kcat values collected from BRENDA and other publications - manually curated
        kcats = pd.read_csv("data/manual_kcat_values.csv")
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
        out_abundance = self.E_data.join(new_abundance, how='left')
        out_abundance.to_csv("cache/enzyme_conc_across_conditions.csv")
        
        new_flux = self._calculate_pFBA(growth_params)
        out_flux = self.V_data.join(new_flux, how='left')
        out_flux.to_csv("cache/pFBA_distribution_across_conditions.csv")

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

    def _calculate_pFBA(self, growth_params):
        
        convert_to_irreversible(self.model)            

        cond = growth_params['carbon']
        uptake = growth_params['growth_rate']
        oxy = growth_params['oxygen']
                
        rxns = dict([(r.id, r) for r in self.model.reactions])
        rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
        try:
            rxns['EX_' + cond + '_e'].lower_bound = -1000 # redefine sole carbon source uptake reaction in mmol/gr/h
        except:
            rxns['EX_glc_e'].lower_bound = -1000
        rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = uptake            
        rxns['EX_o2_e'].lower_bound = oxy
        print "solving model for %s..." %growth_params['title'],
        optimize_minimal_flux(self.model, already_irreversible=True)
        print "\t GR: %.02f 1/h" %self.model.solution.f
        
        flux_dist = pd.DataFrame(self.model.solution.x_dict.items()).set_index(0)
        
        #convert units from mmol/gCDW/h to molecules/gCDW/s
        flux_dist = flux_dist * 6.02214129e23 / 1000 / 3600 

        return flux_dist    
                
    def _calculate_pFVA(self, growth_params):
        '''
            calculates the uncertainties of flux by running FVA 
            given minimization of the total sum of fluxes
            
            Arguments:
                growth_params - 
                
            Returns:
                flux_ranges - the range of each flux
        '''
        
        cond = growth_params['carbon']
        uptake = growth_params['growth_rate']
        oxy = growth_params['oxygen']
        
        convert_to_irreversible(self.model)

        rxns = dict([(r.id, r) for r in self.model.reactions])

        rxns['EX_glc_e'].lower_bound = 0 # uptake of carbon source reaction is initialized    
        rxns['EX_' + cond + '_e'].lower_bound = -1000 # redefine sole carbon source uptake reaction in mmol/gr/h

        rxns['EX_o2_e'].lower_bound = oxy

        rxns['Ec_biomass_iJO1366_core_53p95M'].upper_bound = uptake            
        rxns['Ec_biomass_iJO1366_core_53p95M'].lower_bound = uptake            

        
        fake = Metabolite(id='fake')
        model.add_metabolites(fake)        
                
        for r in self.model.reactions:
            r.add_metabolites({fake:1})
            
        flux_counter = Reaction(name='flux_counter')
        flux_counter.add_metabolites(metabolites={fake:-1})                

        model.add_reaction(flux_counter) 
        self.model.change_objective(flux_counter)
        
        print "solving pFVA"
        fva = flux_variability_analysis(self.model, 
                                        reaction_list=self.model.reactions, 
                                        objective_sense='minimize',
                                        fraction_of_optimum=1.1)
            
        flux_ranges = pd.Series(index=map(lambda x: x.id, self.model.reactions))

        for r, v in fva.iteritems():
            flux_ranges[r] = v['maximum'] - v['minimum']

#        flux_ranges.to_csv("cache/reactions_to_pfva_ranges")

        return flux_ranges

    def manually_remove_reactions(self, reactions):
        # PPKr_reverse reaction is used for ATP generation from ADP 
        # in the FBA model. Nevertheless, acording to EcoCyc, it is used to 
        # to generate polyP (inorganic phosphate) chains from ATP and it is not
        # part of the oxidative phosphorilation
        list_of_reactions = ['PPKr_reverse']
        return reactions.drop(list_of_reactions) 
        
if __name__ == "__main__":

    minimal_conditions = 5
    model_fname = "data/iJO1366_curated.xml"
    model = create_cobra_model_from_sbml_file(model_fname)
    rate = RCAT(model)
    growth_params = {'title':'glucose_heinamnn', 'carbon':'glc',
                     'growth_rate':0.6, 'oxygen':-1000}

    kcat = rate.get_kcat_of_model_reactions()
    rcat_max = rate.get_rcat_max(minimal_conditions)
    rcat_second_max = rate.get_rcat_second_max(minimal_conditions-1)
    r = kcat.index & rcat_max.index
    
    r = rate.manually_remove_reactions(r)


    
    from scipy import stats
    import matplotlib.pyplot as plt
    cor, pval = stats.pearsonr(np.log10(kcat[r]), np.log10(rcat_max[r])) 
    print cor**2, 
    cor, pval = stats.pearsonr(np.log10(kcat[r]), np.log10(rcat_second_max[r])) 
    print cor**2
    
    res =  np.log10(kcat[r]) - np.log10(rcat_max[r])
    res.sort()
    
    fig = plt.figure(figsize=(5,5))
    plt.scatter(kcat[r], rcat_max[r], color='r')
    plt.scatter(kcat[r], rcat_second_max[r], color='0.5')
    plt.xscale('log')    
    plt.yscale('log')    
    plt.xlim(1e-4, 1e4)
    plt.ylim(1e-4, 1e4)    


    #pfva_ranges = rate._calculate_pFVA(growth_params)
    

#    def reactions_proteomics(self):
            
#        # read vilu proteomics file
#        vilu = pd.read_csv(vilu_proteomics) 
#        vilu.set_index("bnumber", inplace=True) 
#        vilu = vilu * 1e12 / 0.3 # convert from copies/fl(cytoplasme) to copies/gCDW
#        vilu.dropna(how='all', inplace=True) # copies/gCDW
#        
#        # read heinmann proteomics file
#        hein = pd.read_csv(heinmann_proteomics)
#        hein.set_index("bnumber", inplace=True)
#        hein.dropna(how='all', inplace=True) # copies/cell        
#        hein = hein[hein>1] # filter enzymes genes with very low expression
#        
#        # convert heinmann units from copies/cell to copies/gCDW

#
#        ab = vilu.join(hein, how='outer') # merge data sets
#        
#        # match bnumbers of genes to reactions in model
#        r_to_single_gene = self.map_model_reaction_to_genes().set_index(1)
#        
#        r_to_ab =  r_to_single_gene.join(ab, how='inner')
#        r_to_ab.set_index(0,inplace=True)
#        r_to_ab = r_to_ab[self.growth_conditions]
#        r_to_ab.to_csv("cache/reactions_to_enzyme_concentrations_across_conditions.csv")
#    
#        return r_to_ab
#
#    def pfba_fluxes(self):
#
#        model = create_cobra_model_from_sbml_file(model_fname)
#        
#        out_flux = pd.DataFrame(columns = self.growth_conditions, index=model.reactions)
#        for cond in self.growth_conditions:
#
#            # load cobra model
#            model = create_cobra_model_from_sbml_file(model_fname)
##        

#_get_pFBA
#    def reactions_dG(self,c = 1e-3, pH=7.5, I=0.2, T=298.15):
#    
#        if not os.path.exists(CC_CACHE_FNAME):
#            cc = ComponentContribution()
#            cc.save_matfile(CC_CACHE_FNAME)
#        else:
#            cc = ComponentContribution.from_matfile(CC_CACHE_FNAME)    
#            
#        reaction_strings, reaction_CID_sparse = self.reactions_formula_and_sparse()
#
#
#        reactions = []
#        reac_strings = []
#        for key, val in reaction_strings.iteritems():    
#            reactions.append(key)
#            reac_strings.append(val)
#    
#        Kmodel = KeggModel.from_formulas(reac_strings)
#        Kmodel.add_thermo(cc)
#        dG0_prime, dG0_std = Kmodel.get_transformed_dG0(pH=7.5, I=0.2, T=298.15)
#        # use 1mM as the standard concentration
#        print c
#        conc = np.ones((1, len(Kmodel.cids))) * 1e-3
#
#        for cid in ['C00001']: # set the conc of water as 1M
#            i = Kmodel.cids.index(cid)
#            conc[0, i] = 1
#            
#        r_to_dGc = pd.DataFrame(index = reactions, columns = ['glc', 'glyc', 'ac', 'std'])
#        for condition in ['glc', 'glyc', 'ac']:
#            for ID, cid  in self.known_cids.iteritems(): # set the conc of known metabolites
#                if ID in Kmodel.cids:
#                    i = Kmodel.cids.index(ID)
#                    c = self.metabolites_concentration[condition][cid]
#                    if not np.isnan(c):
#                        conc[0, i] = c
#                
#            dGc_prime = dG0_prime + R * default_T * np.dot(np.log(conc), Kmodel.S).T
#            r_to_dGc[condition] = dGc_prime
#
#        r_to_dGc['std'] = dG0_std
#        r_to_dGc.to_csv('cache/reactions_to_dGc.csv')
#
#
#class CATALYTIC_RATES(CACHABLE):
#
#    # the minimal set of conditions in which an enzyme is active to allow the calculation of the maximal catalytic rate in vivo
#    minimal_conditions = 5
#
#    def reactions_fluxes(self):
#        pFBA = pd.read_csv('cache/reactions_to_pFBA_across_conditions.csv').set_index('Unnamed: 0')
#        return pFBA * 6.02214129e23 / 1000 / 3600 #convert units from mmol/gCDW/h to molecules/gCDW/s
#        
#    def reactions_abundances(self):
#        abundances = pd.read_csv('cache/reactions_to_enzyme_concentrations_across_conditions.csv').set_index('0')
#        return abundances
#        

#

#        
#    def BRENDA_kcat(self):
#
#        fi = csv.DictReader(open(kcat_fname))
#        ec_to_r = {}
#        for row in fi:
#            if row['Kcat(1/s)'] != 'NaN':
#                r = row['reaction name'].strip()
#                ec_to_r[row['EC number']] = r
#
#        csv_reader = csv.DictReader(open(brenda, 'r'))
#        ec_to_kcat = {}
#        for row in csv_reader:
#            if row["Organism ID"] == '6':
#                if row["kcat"] != "NaN":
#                    kcat = float(row["kcat"])
#                    ec = ".".join([row["EC1"],row["EC2"],row["EC3"],row["EC4"]])
#                    if ec not in ec_to_kcat:
#                        ec_to_kcat[ec] = [kcat]        
#                    else:
#                        ec_to_kcat[ec].append(kcat)
#
#
#        r_to_kcat = {ec_to_r[ec]:gmean(np.array(kcat)) for ec,kcat in ec_to_kcat.iteritems() 
#                    if ec in ec_to_r}
#        return pd.Series(r_to_kcat, index=r_to_kcat.keys())
#        
#    def maximal_rates(self):
#        
##        bnumber = pd.Series(self.map_model_reaction_to_genes().set_index(0))
##        print bnumber
#        brenda = self.BRENDA_kcat()
#        kcats = self.reactions_to_kcat()
#
#        rates = self.enzyme_rates()
#        
#        # return the maximal rate of enzymes and the best condition for each enzyme        
#        max_rates = pd.concat([brenda, kcats, rates.max(axis=1),(rates.idxmax(axis=1))], axis=1)
#        max_rates.columns = ['BRENDA_kcat', 'kcat_vitro', 'rcat_max', 'best_condition']
#        return max_rates 
#
#    def flux_uncertainty(self):
#        v = pd.read_csv('cache/reactions_to_pFBA_across_conditions.csv').set_index('Unnamed: 0').astype('float')
#        fva = pd.read_csv("cache/reactions_to_fva_range.csv").set_index('Unnamed: 0').astype('float')
#        
#        rcatmax = self.maximal_rates()
#        cmax = rcatmax['best_condition'].dropna()
#
#        uncertainty = (fva / v).replace([np.inf, -np.inf], np.nan)
#        
#        out = pd.Series()
#        for r in cmax.index:
#            out[r] = uncertainty[cmax[r]][r]
#        return out
    

        

#class MM_KINETICS(CATALYTIC_RATES):
#    x = 5
#        
#    def reactions_formula_and_sparse(self):
#
#        FORMULA_CACHE = 'cache/reactions_formula.json'
#        REACTION_CID_SPARSE = 'cache/reactions_CID_sparse.json'
#        
#        if not os.path.exists(FORMULA_CACHE) or not os.path.exists(REACTION_CID_SPARSE):
#            reaction_strings, reaction_CID_sparse = self.generate_reactions_formula_sparse()
#        else:            
#            with open(FORMULA_CACHE) as fp:
#                    reaction_strings = json.load(fp)
#    
#            with open(REACTION_CID_SPARSE) as fp:
#                    reaction_CID_sparse = json.load(fp)
#
#        return reaction_strings, reaction_CID_sparse
#
#    def reactions_to_EC(self):
#
#        REACTIONS_EC = 'cache/reactions_to_ec.json'
#        
#        if not os.path.exists(REACTIONS_EC):
#            reactions_EC = self.generate_reactions_to_EC_cache()
#        else:            
#            with open(REACTIONS_EC) as fp:
#                    reactions_EC = json.load(fp)
#        return reactions_EC
#
#    def substrates_km(self):
#    
#        r_strings, r_CID_sparse = self.reactions_formula_and_sparse()
#
#        for r,v in r_CID_sparse.iteritems():
#            r_CID_sparse[r] = {int(cid[1:]):i for cid,i in v.iteritems()}        
#
#        # map reactions to substrates
#        # substrates have positive stoichometric coefficient in the sparse
#        r_to_substrates = {}
#        for r, v in r_CID_sparse.iteritems():
#            r_to_substrates[r] = {i:j for i,j in v.iteritems() if j<0} 
#        
#        r_to_ec = self.reactions_to_EC()
#        ec_to_sparse = {r_to_ec[r]:v for r,v in r_CID_sparse.iteritems()}
#    
#        br = csv.DictReader(open(brenda, 'r'))
#        ec_to_km_sub = defaultdict(dict)
#        for row in br:
#            if row['KM']!='NaN' and row['Organism ID'] == '6':
#                ec = '.'.join([row['EC1'], row['EC2'], row['EC3'], row['EC4']])
#                if ec in ec_to_sparse:
#                    cid = int(row['Compound Id (KEGG)'])
#                    if cid in ec_to_sparse[ec]:
#                        direction = ec_to_sparse[ec][cid]
#                        if direction < 0 :#compound is substrate
#                            ec_to_km_sub[ec][cid] = float(row['KM']) / 1000.0 # convert km units from uM to 
#                
#        r_to_kms = defaultdict(dict)
#        for r,v in r_to_substrates.iteritems():
#            for ec, kms in ec_to_km_sub.iteritems():
#                if r_to_ec[r] == ec:
#                    for c,s in kms.iteritems():
#                        if c in v:
#                           r_to_kms[r][c] = s
#    
#        return {r:v for r,v in r_to_kms.iteritems() if len(v) == len(r_to_substrates[r])}
#

#
#    def backwards_reaction_effect(self):
#
#        dG = pd.read_csv("cache/dGc_for_model_reactions.csv")        
#        dG = dG.set_index(['Unnamed: 0'])
#        
#        thermo = dG[['glc', 'glyc', 'ac']]
#        thermo = 1-np.exp(thermo/(R*default_T))
#        thermo = thermo.replace([np.inf, -np.inf], np.nan)
#
#        return dG, thermo
#
#    def under_saturation_effect(self, include_cofactors=True, include_metabolomics=True):
#
#        reaction_strings, reaction_CID_sparse = self.reactions_formula_and_sparse()
#
#        r_to_saturation = pd.DataFrame(index = self.substrates_km().keys(), columns = ['glc', 'glyc', 'ac'])
#        for condition in ['glyc', 'ac', 'glc']:
#            unknown_reactions = set([])
#            for r, sparse in self.substrates_km().iteritems():
#                si = 1
#                for cid, km in sparse.iteritems():
#                    coef = -1*reaction_CID_sparse[r]['C'+'%05i'%cid]
#                    if cid in self.known_cids.values():
#                        conc = self.metabolites_concentration[condition][cid]
#                        if np.isnan(conc):
#                            unknown_reactions.add(r)
#                            conc = 0.1
#                    else:
#                        unknown_reactions.add(r)
#                        conc = 0.1
##                    s = (km/conc)**coef
#                    s = (conc/km)**coef
#                    si *= s
#                r_to_saturation[condition][r] = si/(1+si)
#        return r_to_saturation, unknown_reactions 
#
#    def damping_factor(self):
#
#        s, unknown_reactions = self.under_saturation_effect()
#        dG, t = self.backwards_reaction_effect()
#        
#        rmaxs, cmaxs = self.maximal_rates()
#        reactions = cmaxs.index & s.index & t.index & self.reactions_to_kcat().index
#    
#        damp = s.loc[reactions]*t.loc[reactions]
#        minimal_damp = pd.Series()
#        max_s = pd.Series()
#        max_t = pd.Series()
#        for reac, cond in cmaxs.loc[reactions].iteritems():
#            if cond in ['glc', 'ac', 'glyc']:
#                minimal_damp[reac] = damp[cond][reac]
#                max_s[reac] = s[cond][reac]
#                max_t[reac] = t[cond][reac]
#            elif cond in ['fum', 'pyr']:
#                minimal_damp[reac] = damp['ac'][reac]
#                max_s[reac] = s['ac'][reac]
#                max_t[reac] = t['ac'][reac]
#            else:
#                minimal_damp[reac] = damp['glc'][reac]
#                max_s[reac] = s['glc'][reac]
#                max_t[reac] = t['glc'][reac]
#        return damp, minimal_damp, max_s, max_t


    
#    model_fname = "data/iJO1366_curated.xml"
#    model = create_cobra_model_from_sbml_file(model_fname)
#    M = MODEL(model, ')
#    r = model.reactions.get_by_id('RPI')
#    sparse = M.map_model_reactions_to_EC(model_fname)
#    upids = {row[0:5]:row[48:54] for row in open(ecoli_genes, 'r')}
#    sys.path.append(os.path.expanduser('~/git/biopython'))
#    from Bio import ExPASy
#    from Bio import SeqIO 
#    from Bio.SeqUtils.ProtParam import ProteinAnalysis
#    
#    C = CACHABLE()
#    reactions = C.map_model_reaction_to_genes().set_index(0)
#    enzyme_data = pd.DataFrame(columns=["bnumber", "gene", "MW [KDa]"],
#                       index=reactions.index)
#    enzyme_data.index.name = "reaction"
#    
#    for i, r in enumerate(reactions.index):
#        b = reactions.loc[r][1]
#
#        enzyme_data["bnumber"][r] = b
#        
#        try:
#            enzyme_data["gene"][r] = C.gene_names[b]  
#            handle = ExPASy.get_sprot_raw(upids[b])
#            my_seq = str(SeqIO.read(handle, "swiss").seq)
#        except:
#            continue
#        print i, r
#
#        enzyme_data["MW [KDa]"][r] = ProteinAnalysis(my_seq).molecular_weight()/1000.0
#    
#    enzyme_data.to_csv("cache/enzymatic_data.csv")
#    CR = CATALYTIC_RATES()
#    fva = CR.pFVA()
#    kcat = CR.maximal_rates()['kcat_vitro']
#    rcatmax = CR.maximal_rates()['rcat_max']
#    
#    x = (rcatmax / kcat).dropna()
#    y = CR.flux_uncertainty()[x.index]
#    
#    import matplotlib.pyplot as plt    
#    fig = plt.figure()
#    ax = plt.axes()
#    plt.plot(x,y, 'ro')
#    
#    ax.set_xscale('log')
#    ax.set_yscale('log')
    
 
#growth_rates = {'ac':0.29, 'fum':0.47, 'gal':0.17, 'glc':0.60, 'gam':0.39, 'glyc':0.47, 'pyr':0.4, 'succ':0.4, 
#                  'anaerobic':0.55, 'Chem_05':0.5, 'Chem_035':0.35, 'Chem_02':0.2, 'Chem_012':0.12,
#                  '42Cdeg':0.65, 'pH6':0.5, 'NaCl_50mM':0.65, 
#                  'vilu_011':0.11, 'vilu_021':0.21, 'vilu_031':0.31, 'vilu_04':0.40, 'vilu_049':0.49}
#
#cond_to_cell_vol = {'ac':2.4e-15, 'fum':2.4e-15, 'gal':1.9e-15, 'glc':3.2e-15, 'gam':2.9e-15, 'glyc':2.3e-15, 
#	               'pyr':2.1e-15, 'succ':2.4e-15, 'anaerobic':2.9e-15, 'Chem_05':2.6e-15, 'Chem_035':2.4e-15, 
#                    'Chem_02':2.2e-15, 'Chem_012':2.1e-15, '42Cdeg':2.8e-15, 'pH6':3.1e-15, 'NaCl_50mM':2.8e-15,
#                    'vilu_011':0.73e-15, 'vilu_021':0.91e-15, 'vilu_031':1.22e-15, 'vilu_04':1.46e-15, 'vilu_049':1.69e-15} # in liter
#
## sort conditions by growth rate
#growth_conditions = sorted(growth_rates.keys(), key=growth_rates.get)    
## imported files
#
#metabolites_concentration = pd.read_csv(metabolites_conc) 
#metabolites_concentration.set_index('Compound Id (KEGG)', inplace=True)
#metabolites_concentration.dropna(how='all', inplace=True)
#known_cids = {'C'+'%05i'%c:c for c in metabolites_concentration.index}                
#    
#model_fname = "data/iJO1366_curated.xml"
#ecoli_genes = "data/all_ecoli_genes.txt"
#heinmann_proteomics = "data/Heinmann_proteomics.csv"
#vilu_proteomics = "data/Vilu_proteomics.csv"
#brenda = "data/brenda_data.csv"
#metabolites_conc = "data/metab_conc.csv"
#    genes = [row[0:5] for row in open(ecoli_genes, 'r')]
#    gene_names = {row[0:5]:row[84:].split(';')[0].strip() for row in open(ecoli_genes, 'r')}

#CC_CACHE_FNAME = os.path.expanduser('~/git/component-contribution/cache/component_contribution.mat')
   
    
    
