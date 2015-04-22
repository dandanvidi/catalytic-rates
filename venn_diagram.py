from kapp import MM_KINETICS, PLOT_DATA
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)
mm = MM_KINETICS(model)

KMs = pd.DataFrame.from_csv('../cache/KM_values.csv')
kcat = mm.get_kcat_of_model_reactions()
rcat = mm.calculate_enzyme_rates()
rmax = mm.get_rcat_max(7)


flux = mm.V_data
expression = mm.E_data
kcats = set(kcat.index)
rmaxs = set(kcat.index)

all_reactions = set(map(lambda x: x.id, model.reactions))

#plt.figure()
#venn3([expression, kcats, rmaxs], set_labels=('reactions', 'kcat', 'rmax'))