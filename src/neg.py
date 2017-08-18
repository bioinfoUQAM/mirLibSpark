'''
script to create negative miRNA in raw format

author: Chao-Jung Wu
date: 2017-08-17
'''
import utils as ut
import random



def collect_negtive_miRNA ():
  known_ath100 = []
  #f_known = '/home/cloudera/Desktop/mirLibHadoop/test/high_conf_mature_ath_uniq.txt'
  f_known = '/home/cloudera/Desktop/mirLibHadoop/test/mirbase_all_427_real_ath.txt'
  with open (f_known, 'r') as fh:
    for i in fh:
      known_ath100.append(i.rstrip('\n').split('\t')[0])

  genome_path = '/home/cloudera/workspace/mirLibHadoop/dbs/ATH/Genome/'
  file_ext = '.fas'
  genome = ut.getGenome (genome_path, file_ext)
  chr1seq = genome['Chr1']
  lenchr1 = len(chr1seq)
  negseqs = []

  i = 1
  while i < 1001:
    start = random.randint(0, lenchr1-25)
    end = random.randint(start+17, start+24)
    length = end - start + 1
    seq = chr1seq[start: end+1]
    if seq[:10] in known_ath100: continue
    if seq[-10:] in known_ath100: continue
    if seq not in negseqs: 
      negseqs.append(seq)
      i += 1
  return negseqs

negseqs = collect_negtive_miRNA ()
outfile = '../negativeSet/neg_ath1000.txt'
fh_out = open (outfile, 'w')
for s in negseqs:
  print >> fh_out, s, '\t', 500
fh_out.close()



