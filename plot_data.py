from kapp import MODEL 


import numpy as np

from cobra.io.sbml import create_cobra_model_from_sbml_file

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
M = MODEL(model)


