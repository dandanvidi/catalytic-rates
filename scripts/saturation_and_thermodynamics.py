import sys, os, csv, ast
import uncertainties.unumpy as unumpy  
from uncertainties import ufloat
from rcat import RCAT
import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np
sys.path.append(os.path.expanduser('~/git/thermodynamics-for-cobra'))
from thermodynamics_for_cobra import THERMODYNAMICS_FOR_COBRA
from component_contribution.kegg_reaction import KeggReaction
from component_contribution.kegg_model import KeggModel
from component_contribution.component_contribution import ComponentContribution
from component_contribution.thermodynamic_constants import R, default_T

class MM_KINETICS(THERMODYNAMICS_FOR_COBRA, RCAT):

    def __init__(self, model, reactions):
        RCAT.__init__(self)

        metabolites = pd.DataFrame.from_csv("../data/model_metabolites_to_kegg.txt", sep='\t')
        THERMODYNAMICS_FOR_COBRA.__init__(self, model, reactions, metabolites)
        
        self.CC_CACHE_FNAME = os.path.expanduser('../../component-contribution/cache/component_contribution.mat')

        self.metab_conc = pd.read_csv('../data/bennett_metabolite_concentrations[mM].csv') # concentrations in mM
        self.metab_conc.set_index('CID', inplace=True) 
        self.metab_conc.dropna(how='all', inplace=True)
        
        self.known_cids = {'C'+'%05i'%c:c for c in self.metab_conc.index}
        
        self.substrates = self.get_reaction_substrates()
        self.Ks = self.get_known_Ks_values()

    def set_metabolite_concentrations(self, condition):        
        #initialize concentrations of metabolites as 100uM
        conc = np.ones((1, len(self.Kmodel.cids))) * 1e-4 # concentrations in M
        for ID, cid  in self.known_cids.iteritems(): # set the conc of known metabolites
            if ID in self.Kmodel.cids:
                i = self.Kmodel.cids.index(ID)
                # set the conc of water as 1M
                if ID == 'C00001':
                    conc[0, i] = 1
                else:
                    c = self.metab_conc[condition][cid]
                    if not np.isnan(c):
                        conc[0, i] = c  * 1e-3  # concentrations in M

        return conc
    
    def get_udGm_prime(self):
        conc = np.ones((1, len(self.Kmodel.cids))) * 1e-3 # concentrations in M
        udGm_prime = self.udG0_prime + R * default_T * np.dot(np.log(conc), self.Kmodel.S)
        udGm_prime[udGm_prime>200] = unumpy.uarray(200, 0)
        udGm_prime[udGm_prime<-200] = unumpy.uarray(-200, 0)        
        return np.array(udGm_prime)[0]
        
    def get_udGc_prime(self, condition):
        conc = self.set_metabolite_concentrations(condition)
        udGc_prime = self.udG0_prime + R * default_T * np.dot(np.log(conc), self.Kmodel.S)
        
        return np.array(udGc_prime)[0]

    def get_thermodynamic_effect(self, condition):
        udGc_prime = self.get_udGc_prime(condition)
        udGc_prime[udGc_prime>0] = 0
        return -unumpy.expm1(udGc_prime/(R*default_T))

    def get_known_Ks_values(self):
        metabolites = self.metabolites
        metabolites.set_index('name', inplace=True)
        cids = metabolites.kegg_id.dropna().astype('str').to_dict()
        km_data = csv.reader(open("../data/KM_values_from_ecocyc.csv", 'r'), delimiter='\t')
        km_data.next()
        km_sparses = {}
        for row in km_data:
            r = row[0][1:-1]
            substrates = ast.literal_eval(row[1][1:-1])
            km_sparse = {}
            for s, ks in substrates.iteritems():
                s = s.replace("*", "'")
                if ks !=[] and s in cids:
                    ks = unumpy.uarray(np.mean(ks), np.std(ks)) / 1000 # in mM
                    km_sparse[cids[s]] = ks
                elif ks !=[] and s not in cids:
                    print s
                km_sparses[r] = km_sparse
        Ks = [km_sparses[k] if k in km_sparses.keys() else {} for k in self.reactions]
        return Ks

    def get_reaction_substrates(self):
        substrates_list = []
        for i, r in enumerate(self.reactions):
            substrates = [k for k, v in self.reaction_sparses[i].iteritems() if v < 0]
            substrates_list.append(substrates)
        return substrates_list
        
    def get_prod_s_over_Ks(self, condition):
        s_over_Ks = unumpy.uarray([1]*len(self.reactions), [0]*len(self.reactions))
        for i, r in enumerate(self.reactions):
            if set(self.substrates[i]).issubset(set(self.Ks[i].keys())):
                for k, v in self.Ks[i].iteritems():
                    if k not in self.known_cids:
                        c = 0
                    else:
                        c = self.metab_conc[condition][self.known_cids[k]]
                    ms = np.abs(self.reaction_sparses[i][k])
                    s_over_Ks[i] *= (c / v) ** ms
            else:
                s_over_Ks[i] = unumpy.uarray(0, 0)
        return s_over_Ks
    
    def get_saturation_effect(self, condition):
        s_over_Ks = self.get_prod_s_over_Ks(condition)
        return s_over_Ks / (1+s_over_Ks)
    
if __name__ == "__main__":
    rate = RCAT()
    mm = MM_KINETICS(rate.model, ['MDH'])

