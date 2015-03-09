import csv
import pandas as pd
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
import ast
from scipy.stats.mstats import gmean
import numpy as np

model_fname = "../data/iJO1366_curated.xml"
model = create_cobra_model_from_sbml_file(model_fname)
convert_to_irreversible(model)

reactions = set([r.id for r in model.reactions])
metabolites = set([m.name for m in model.metabolites])

kms = csv.reader(open("../data/Km_values.csv",'r'), delimiter='\t')
kms.next()
data = pd.DataFrame(columns=metabolites, index=reactions)
for i, row in enumerate(kms):
    r, s , p = row[0][1:-1], row[1][1:-1], row[2][1:-1]
    if r == 'PGI':
        print s
    s = ast.literal_eval(s)
    for m, ks in s.iteritems():
        if ks != []:
            m = m.replace("*", "'")
            data.loc[r,m] = -gmean(np.array(ks))

    p = ast.literal_eval(p)
    for m, kp in p.iteritems():
        if kp != []:
            m = m.replace("*", "'")
            data.loc[r,m] = gmean(np.array(kp))

data.dropna(how='all', inplace=True)
data.to_csv('../cache/KM_values.csv')