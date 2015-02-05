''' saturation and thermodynamics '''

from kapp import MM_KINETICS
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)
mm = MM_KINETICS(model)

kcat = mm.get_kcat_of_model_reactions()
rcat = mm.calculate_enzyme_rates()
rmax = mm.get_rcat_max(7)

mm.cache_reactions_dG()
#dG, T = mm.backwards_reaction_effect()
#
#reactions = rmax.index & kcat.index & T.index
#kcat = kcat[reactions]
#rmax = rmax[reactions]
#T = T.loc[reactions]
#
#growth_rates = rc.growth_conditions.growth_rate_1_h
#
#
#damp, minimal_damp, s, t = damping_factor()
#reactions = minimal_damp.index
##dist = residual(np.log10(kcats[r]), np.log10(kmaxs[r]))[reactions]
#dist = np.log10(kmaxs/kcats)
#S = SATURATION()
#sat, unknown_reactions = S.saturation_index()
#T = THERMODYNAMICS() 
#dG, ter = T.reversibility_index()
#unknown = unknown_reactions & set(reactions)
#known = set(reactions) - unknown