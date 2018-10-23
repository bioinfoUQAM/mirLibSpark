'''
script to create negative miRNA in raw format

author: Chao-Jung Wu
date: 2017-08-17
'''
import utils as ut
import random

def collect_negtive_miRNA ():
  known_ath100 = []
  f_known = '../input_samples/mirbase_all_349_real_ath_uniq.txt'
  with open (f_known, 'r') as fh:
    for i in fh:
      known_ath100.append(i.rstrip('\n').split('\t')[0])

  genome_path = '../dbs/ATH_TAIR10/Genome/'
  file_ext = '.fa'
  genome = ut.getGenome (genome_path, file_ext)
  genomeseq_all = ''
  for k, v in genome.items():
    genomeseq_all += v + ut.getRevComp (v)
  lengenome = len(genomeseq_all)
  negseqs = []

  i = 1
  while i < 1001:
    start = random.randint(0, lengenome-25)
    end = random.randint(start+17, start+24)
    length = end - start + 1
    seq = genomeseq_all[start: end+1]
    if seq[:10] in known_ath100: continue
    if seq[-10:] in known_ath100: continue
    if seq not in negseqs: 
      negseqs.append(seq)
      i += 1
  return negseqs

negseqs = collect_negtive_miRNA ()
outfile = 'neg_ath1000.txt'
fh_out = open (outfile, 'w')
for s in negseqs:
  print >> fh_out, s+'\t500'
fh_out.close()



