import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np

class RCAT(object):
    
    def __init__(self):
        
        self.path = '/home/dan/Dropbox/PhD/catalytic_rates/Figures'
        self.model_fname = "../data/iJO1366.xml"
        self.model = create_cobra_model_from_sbml_file(self.model_fname)
        convert_to_irreversible(self.model)  
        self.rxns = dict([(r.id, r) for r in self.model.reactions])
        growth_conditions = "../data/growth_conditions.csv"
        self.gc = pd.DataFrame.from_csv(growth_conditions)

        self.v = pd.DataFrame.from_csv('../cache/flux[mmol_h_gCDW].csv')
        self.v.replace(0, np.nan, inplace=True)
        self.v.dropna(how='all', inplace=True)
        
        self.p = pd.DataFrame.from_csv("../cache/abundance[copies_fl].csv")
        self.p = self.p[self.p>=10]# remove genes below 10 copies/fL
        self.p.dropna(how='all', inplace=True)
        
        self.enzymatic_reactions = self.enzymatic_reactions()       
        
        self.kcat = self.get_kcat()
        self.rcat = self.calculate_catalytic_rates()
        self.v_over_E = self.get_v_over_E()        
        self.minimal_conditions = 5
        self.rmax = self.get_rmax()
        self.rcatn = self.get_rcatn()
        self.rmaxn = self.get_rmaxn()
        self.vmax = self.get_Vmax() 
        
    def map_reactions_to_gene_names(self):
        gene_names = {row[0:5]:row[84:].split(';')[0].strip() 
                      for row in open("../data/all_ecoli_genes.txt", 'r')}
        return {k:gene_names[v] for k,v in 
                self.reactions_to_homomeric_enzymes().iteritems()
                if v in gene_names.keys()}
            
    def convert_copies_fl_to_copies_gCDW(self, proteomics):
        rho = 1100 # average cell density gr/liter
        DW_fraction = 0.3 # fraction of DW of cells
        proteomics = proteomics * 1e15 / (rho * DW_fraction)
        return proteomics

    def convert_copies_gCDW_to_copies_fl(self, proteomics):
        rho = 1100 # average cell density gr/liter
        DW_fraction = 0.3 # fraction of DW of cells
        proteomics = proteomics * (rho * DW_fraction) / 1e15
        return proteomics

    def convert_mmol_gCDW_h_to_molecules_gCDW_s(self, flux):
        flux = flux * 6.02214129e23 / 1000/ 3600  
        return flux
    
    def convert_kDa_fl_to_mg_gCDW(self, array):
        rho = 1100 # average cell density gr/liter
        DW_fraction = 0.3 # fraction of DW of cells
        return array * 1.66053892e-3 / (rho * DW_fraction)
        
    def enzymatic_reactions(self):
        return filter(lambda r:len(r.genes)>=1, self.model.reactions)
        
    def reactions_to_unique_enzyme(self):

        one_enzyme_reac = filter(lambda r: 'or' not in r.gene_reaction_rule, 
                                 self.enzymatic_reactions)
        one_enzyme_reac = {r.id:map(lambda g: g.id, r.genes) for r in one_enzyme_reac}   
        
        one_enzyme_reac = self.include_known_isoenzmyes(one_enzyme_reac)
        return one_enzyme_reac
        
    def include_known_isoenzmyes(self, dictionary):
        dictionary['METS'] = ["b3829"] # metE - cobalamin-independent homocysteine transmethylase
        dictionary["HCO3E"] = ["b0126"] # can - carbonic anhydrase
        dictionary["PFK"] = ["b3916"] #  6-phosphofructokinase                
        dictionary["RPI"] = ["b2914"] #  ribose-5-phosphate isomerase A
        return dictionary
        

    def reactions_to_enzyme_complex(self):
        return {k:v for k,v in self.reactions_to_unique_enzyme().iteritems()
                                if len(v)>1}


    def reactions_to_homomeric_enzymes(self):
        ''' 
            find and match all reactions in the model 
            that are catalyzed by a single gene 
        '''
        return {k:v[0] for k,v in self.reactions_to_unique_enzyme().iteritems()
                                if len(v)==1}

    def calculate_catalytic_rates(self):
        '''
            calculates the catalytic rate of an enzyme for a given reaction
            by dividing the flux through the reaction by the enzyme copy number. 
            
            Important: 
                Data must have matching units
        '''
        
        p = self.convert_copies_fl_to_copies_gCDW(self.p)
        v = self.convert_mmol_gCDW_h_to_molecules_gCDW_s(self.v)   
        
        rcat = pd.DataFrame(index=self.reactions_to_homomeric_enzymes().keys(), 
                            columns=self.gc.index)
        for reaction, gene in self.reactions_to_homomeric_enzymes().iteritems():
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

        self.rcat.dropna(thresh=self.minimal_conditions, inplace=True)

        rmax = pd.DataFrame(index=self.rcat.index)
        reactions = map(self.model.reactions.get_by_id, rmax.index)
        subsystems = map(lambda r: r.subsystem, reactions)
        
        rmax['rmax [s^-1]'] = self.rcat.max(axis=1)
        rmax['subsystem'] = subsystems
        rmax['condition'] = self.rcat.idxmax(axis=1)
        rmax['carbon source'] = map(lambda x: self.gc['carbon source'][x], rmax.condition)
        rmax.to_csv('../cache/rmax_values.csv')
        return rmax
        
    def get_v_over_E(self):
        
        genes = pd.DataFrame.from_csv('../data/model_genes_to_mass.csv', sep='\t')
        mass = genes['mass(Da)'][self.p.index]/1000 #mass in kDa
        weighted_mass = self.p.multiply(mass, axis=0)
        weighted_mass.dropna(how='all', inplace=True)
        weighted_mass = self.convert_kDa_fl_to_mg_gCDW(weighted_mass)

        reactions = map(lambda x: x.id, self.enzymatic_reactions)
        v_over_E = pd.DataFrame(index=reactions, columns=self.gc.index)

        for r in self.enzymatic_reactions:
            genes = map(lambda x: x.id, r.genes)
            try:
                v_over_E.loc[r.id] = self.v.loc[r.id] / weighted_mass.loc[genes].sum()
            except KeyError:
                continue
        
        v_over_E.replace([0, np.inf, -np.inf], np.nan, inplace=True)
        v_over_E.dropna(how='all', inplace=True)
        v_over_E.drop('PPKr_reverse', axis=0, inplace=True)
        
        return v_over_E * 1000 / 60
        
    def get_Vmax(self):

        self.v_over_E.dropna(thresh=self.minimal_conditions, inplace=True)

        Vmax = pd.DataFrame(index=self.v_over_E.index)
        reactions = map(self.model.reactions.get_by_id, Vmax.index)
        subsystems = map(lambda r: r.subsystem, reactions)

        Vmax['Vmax [umol/mg/min]'] = self.v_over_E.max(axis=1)        
        Vmax['subsystem'] = subsystems
        Vmax['condition'] = self.v_over_E.idxmax(axis=1)
        Vmax['carbon source'] = map(lambda x: self.gc['carbon source'][x], Vmax.condition)
        Vmax.to_csv('../cache/Vmax_values.csv')
        return Vmax

    def get_rcatn(self):
        reactions = self.kcat.index & self.rmax.index
        su_per_as = (self.kcat['subunits'] / self.kcat['Catalytic Sites'])[reactions]
        return self.rcat.loc[reactions].mul(su_per_as, axis=0)

    def get_rmaxn(self):

        su_per_as = self.kcat['subunits'] / self.kcat['Catalytic Sites']

        rmaxn = self.rmax
        rmaxn['rmax [s^-1]'] = self.rmax['rmax [s^-1]'][self.kcat.index]*su_per_as
        return rmaxn[rmaxn['rmax [s^-1]']>0]

    def get_kcat(self):
        
        # kcat values collected from BRENDA and other publications - manually curated
        kcats = pd.read_csv("../data/curated_kcat_values.csv")
        kcats.set_index('reactions', inplace=True)

        return kcats        

if __name__ == "__main__":
    R = RCAT()
