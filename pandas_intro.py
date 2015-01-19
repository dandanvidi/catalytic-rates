# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 12:51:27 2015

@author: access

practice Pandas Data structure

"""
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import numpy as np
import pandas as pd

model_fname = "data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)

conditions = pd.read_csv('data/growth_conditions.csv').title
conditions.name = 'condition'
index1 = [r.id for r in model.reactions for i in (0, 1)]
index2 = []
for i in range(len(model.reactions)+1):
    index2 = index2+['min', 'max']

index = [np.array(index1), np.array(index2)]

data = pd.DataFrame(index=index, columns=conditions)

data.to_csv('cache/pFVA_ranges_template.csv')


