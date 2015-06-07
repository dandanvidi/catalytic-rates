from rcat import RCAT
from saturation_and_thermodynamics import MM_KINETICS
import pandas as pd

R = RCAT()
kcat = R.kcat
rcat = R.rcat
mm = MM_KINETICS(R.model, list(R.rmax.index))
kcat = kcat[mm.reactions]
rmax = R.rmax.loc[mm.reactions]
rmax['carbon source'] = map(lambda x: R.gc['carbon source'][x], rmax.condition)

conditions = ['glc', 'ac', 'glyc']
s = pd.DataFrame(columns=conditions, index=mm.reactions)
t = pd.DataFrame(columns=conditions, index=mm.reactions)

for i, c in enumerate(conditions):
    s[c] = mm.get_saturation_effect(c)
    t[c] = mm.get_thermodynamic_effect(c)
    
s.to_csv('../cache/saturation_effect.csv')
t.to_csv('../cache/thermodynamic_potential.csv')