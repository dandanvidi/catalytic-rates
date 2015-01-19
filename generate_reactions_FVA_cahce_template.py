# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 12:51:27 2015

@author: access

generate pFVA cahce template

"""
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np
import pandas as pd

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)

conditions = list(pd.read_csv('data/growth_conditions.csv').title.values)

index = [(r.id, m) for r in model.reactions for m in ('minimum', 'maximum')]
index += [('flux_counter', 'minimum'), ('flux_counter', 'maximum')]
data = pd.DataFrame(index=index, columns=conditions)

data.to_csv('cache/pFVA_ranges_template.csv')


