from bs4 import BeautifulSoup
import pandas as pd
from scripts.kapp import CACHABLE
import urllib 

C = CACHABLE()
reactions = C.map_model_reaction_to_genes().set_index(0)

genes = {row[0:5]:row[56:63] for row in open('data/all_ecoli_genes.txt', 'r')
                                        if row[0:5] in reactions.values}
new_dict = {}
for j,(b, EG) in enumerate(genes.iteritems()):                                      
    sock = urllib.urlopen("http://ecocyc.org/ECOLI/NEW-IMAGE?type=GENE&object=%s" %EG) 
    html = sock.read() 
    doc = BeautifulSoup(html)
    classes = doc.findAll('p')                           
    subunits = 1    
    for item in classes:
        title = item.contents[0].strip()
        if title == 'Subunit composition of':
            for s in item.findAll('sub'):
                try:
                    subunits = int(s.contents[0].strip())
                except ValueError:
                    continue
            break

    print j, b, "->", subunits, " subunits"      
    new_dict[b] = subunits

subunits = pd.DataFrame(new_dict.items())
subunits.to_csv("cache/subunits.csv")







#
#
#
#
#
#
#
#
#
#
#
#
#    
#    m = 0
#    try:
#        a = doc.findAll('p')[4]
#    except:
#        continue
#    if 'Subunit composition of' in str(a):
#        try:
#            a = doc.findAll('sub')
#        except:
#            continue
#
#    print a
#    if j >2:
#        break
#        if 'Subunit composition of' in str(a):
#            m = int(str(a.sub).split('<sub>')[1][0])
#            break
#    if m == 0:
#        m =1
#    new_dict[b] = m
#    print j, EG, "->", m, " subunits"    


subunits = pd.DataFrame(new_dict.items())
subunits.to_csv("cache/subunits.csv")