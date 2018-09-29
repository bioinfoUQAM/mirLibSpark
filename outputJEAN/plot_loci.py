
def getRevComp (seq):
  intab = "ACGT"
  outab = "TGCA"
  #= v1 ================================
  from string import maketrans
  trantab = maketrans(intab, outab)
  #= v2 ================================
  #trantab = str.maketrans(intab, outab)
  #=====================================
  n_seq = seq.translate(trantab)
  return n_seq[::-1]
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

####################################################################################################


genome_seq_file = '../dbs/JEAN/Genome/XLOC_126981_6D_REVCOMP.fasta'
sRNAmappings_file = 'demo_loci.txt'

with open (genome_seq_file, 'r') as fh: genseq = [x.rstrip('\n').split('\t') for x in fh.readlines()][1][0]



with open (sRNAmappings_file, 'r') as fh: mappings = [x.rstrip('\n').split('\t') for x in fh.readlines()]

for m in mappings: 
  srna = m[0]
  freq = m[1]
  strd = str(m[3])
  if strd == '-': srna = getRevComp (srna)
  if strd == '+': strd = 'p'
  else: strd = 'n'
  posch= m[5]
  line = str(posch).zfill(4) + '.'*int(posch) + srna + '_' + freq + '_' + strd + '_' + str(len(srna))
  print(line)
print('....' + genseq)

