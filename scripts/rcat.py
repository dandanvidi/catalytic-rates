import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np

class RCAT(object):
    
    def __init__(self):
        
        self.model_fname = "../data/iJO1366.xml"
        self.model = create_cobra_model_from_sbml_file(self.model_fname)
        convert_to_irreversible(self.model)  

        growth_conditions = "../data/growth_conditions.csv"
        self.gc = pd.DataFrame.from_csv(growth_conditions)

        self.v = pd.DataFrame.from_csv('../cache/flux[mmol_h_gCDW].csv')
        self.v.replace(0, np.nan, inplace=True)
        self.v.dropna(how='all', inplace=True)
        
        self.p = pd.DataFrame.from_csv("../cache/abundance[copies_fl].csv")
        self.p.dropna(how='all', inplace=True)
        
        self.rcat = self.calculate_catalytic_rates()
        self.rmax = self.get_rmax()
        self.kcat = self.get_kcat()
        self.rxns = dict([(r.id, r) for r in self.model.reactions])
        
    def map_reactions_to_gene_names(self):
        gene_names = {row[0:5]:row[84:].split(';')[0].strip() 
                      for row in open("../data/all_ecoli_genes.txt", 'r')}
        return {k:gene_names[v] for k,v in 
                self.reactions_by_homomeric_enzymes().iteritems()
                if v in gene_names.keys()}
            
    def convert_copies_fl_to_copies_gCDW(self, proteomics):
        rho = 1100 # average cell density gr/liter
        DW_fraction = 0.3 # fraction of DW of cells
        proteomics = proteomics * 1e15 / (rho * DW_fraction)
        return proteomics
        
    def convert_mmol_gCDW_h_to_molecules_gCDW_s(self, flux):
        flux = flux * 6.02214129e18 / 36  
        return flux
        
    def reactions_by_homomeric_enzymes(self):
        ''' 
            find and match all reactions in the model 
            that are catalyzed by a single gene 
        '''
        one_enzyme_reac = filter(lambda r: len(r.genes)==1, self.model.reactions)
        reacs = map(lambda r: r.id, one_enzyme_reac)
        genes = map(lambda r: list(r.genes)[0].id, one_enzyme_reac)
        r_to_single_gene = dict(zip(reacs, genes))
        # manual additions of genes with only one active isoenzyme 
        r_to_single_gene['METS'] = "b3829" # metE - cobalamin-independent homocysteine transmethylase
        r_to_single_gene["HCO3E"] = "b0126" # can - carbonic anhydrase
        r_to_single_gene["PFK"] = "b3916" #  6-phosphofructokinase                
        r_to_single_gene["RPI"] = "b2914" #  ribose-5-phosphate isomerase A 
        return r_to_single_gene

    def calculate_catalytic_rates(self):
        '''
            calculates the catalytic rate of an enzyme for a given reaction
            by dividing the flux through the reaction by the enzyme copy number. 
            
            Important: 
                Data must have matching units
        '''
        p = self.p[self.p>=10]# remove genes below 10 copies/fl

        p = self.convert_copies_fl_to_copies_gCDW(p)
        v = self.convert_mmol_gCDW_h_to_molecules_gCDW_s(self.v)   
        
        rcat = pd.DataFrame(index=self.reactions_by_homomeric_enzymes().keys(), 
                            columns=self.gc.index)
        for reaction, gene in self.reactions_by_homomeric_enzymes().iteritems():
            try:
                rcat.loc[reaction] = v.loc[reaction] / p.loc[gene]
            except KeyError:
                continue
            
        rcat.replace([0, np.inf, -np.inf], np.nan, inplace=True)
        rcat.dropna(how='all', inplace=True)
        
        # PPKr_reverse reaction is used for ATP generation from ADP 
        # in the FBA model. Nevertheless, acording to EcoCyc, it is used to 
        # to generate polyP (inorganic phosphate) chains from ATP and it is not
        # part of the oxidative phosphorilation
        rcat.drop('PPKr_reverse', axis=0, inplace=True)
        rcat.to_csv('../cache/rcat_values.csv')
        return rcat 

    def get_rmax(self):
        minimal_conditions = 5
        self.rcat.dropna(thresh=minimal_conditions, inplace=True)
        rmax = pd.DataFrame(index=self.rcat.index, columns=['rmax [s^-1]'])
        rmax['rmax [s^-1]'] = self.rcat.max(axis=1)
        rmax['condition'] = self.rcat.idxmax(axis=1)
        rmax.to_csv('../cache/rmax_values.csv')
        return rmax
        
    def get_kcat(self):
        
        # kcat values collected from BRENDA and other publications - manually curated
        kcats = pd.read_csv("../data/curated_kcat_values.csv")
        kcats.set_index('reactions', inplace=True)

        return kcats['kcat per subunit']        

if __name__ == "__main__":
    R = RCAT()