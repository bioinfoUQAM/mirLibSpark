'''
script to create negative miRNA

author: Chao-Jung Wu
date: 2017-08-17
'''
import utils as ut
import random

genome_path = '/home/cloudera/workspace/mirLibHadoop/dbs/ATH/Genome/'
file_ext = '.fas'
genome = ut.getGenome (genome_path, file_ext)
#for k in genome.keys(): print k
chr1seq = genome['Chr1']
lenchr1 = len(chr1seq)
negseqs = []

i = 1
while i < 101:
  start = random.randint(0, lenchr1-25)
  end = random.randint(start+17, start+24)
  length = end - start + 1
  seq = chr1seq[start: end+1]
  if seq not in negseqs: 
    negseqs.append(seq)
    i += 1

print len (negseqs)

