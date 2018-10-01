'''
Chao-Jung Wu
2018-09-26

Goal to collect genome-wide all KEGG pathways and their names of an organism.
Targeting organisms are ath, tae, and two other plants.

'''
import sys, os
from bioservices import KEGG

if __name__ == '__main__' :
  if not len(sys.argv) == 2:
    sys.stderr.write('please provide organim 3-letter code. \nusage: bioservicesKEGG [organim 3-letter code]\nexample: python bioservicesKEGG.py ath')
    sys.exit()

  organism = sys.argv[1]
  rep = '../output/' + organism + '/'
  if not os.path.exists(rep): os.makedirs(rep)
  k = KEGG(verbose=False)
  k.organism = organism
  paths = k.pathwayIds
  for p in paths: 
    path = p.split(':')[1]
    outfile = '../output/' + organism + '/' + path + '.txt'
    with open (outfile, 'w') as fh:
      print(k.get(path), file=fh, flush=True)
