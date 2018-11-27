'''
script to create negative miRNA in raw format

random position in the genome with length ranging from 18 to 25.

author: Chao-Jung Wu
date: 2017-08-17
'''
import utils as ut
import random

def get_nonMirna_list (genome_path):
  infile = '../dbs/ATH_TAIR10/TAIR10_ncRNA_CDS.gff'
  genome = ut.getGenome (genome_path, file_ext, 'All') #= genome[chr] = sequence
  l_non_miRNA = [] #= ['TGGATTTATGAAAGACGAACAACTGCGAAA']
  with open (infile, 'r') as fh:
    for i in fh:
      data = i.split('\t') #= ['Chr1', 'TAIR10', 'CDS', '3760', '3913', '.', '+', '0', 'Parent=AT1G01010.1,AT1G01010.1-Protein;\n']
      chromo   = data[0]
      if chromo == 'ChrC': chromo = 'chloroplast'
      if chromo == 'ChrM': chromo = 'mitochondria'
      begin    = int(data[3])
      end      = int(data[4])
      strand   = data[6]
      seq = genome[chromo][begin:end+1]
      if strand == '-':
        seq = ut.getRevComp (seq)
      l_non_miRNA.append(seq)
  return l_non_miRNA

def get_known_miRNA ():
  known_ath349 = []
  f_known = '../input_samples/mirbase_all_349_real_ath_uniq.txt'
  with open (f_known, 'r') as fh:
    for i in fh:
      known_ath349.append(i.rstrip('\n').split('\t')[0])
  return known_ath349

def collect_negtive_miRNA ():
  genome_path = '../dbs/ATH_TAIR10/Genome/'

  known_ath349 = get_known_miRNA () #
  l_non_miRNA = get_nonMirna_list (genome_path) #

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
    if seq[:10] in known_ath349: continue
    if seq[-10:] in known_ath349: continue
    if seq[:10] in l_non_miRNA: continue
    if seq[-10:] in l_non_miRNA: continue
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



