
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
sRNAmappings_file = 'demo_loci.txt' # bowtie result: seq	freq	nbloci	strand	chromo	poschr
with open (genome_seq_file, 'r') as fh: 
  genseq = [x.rstrip('\n').split('\t') for x in fh.readlines()][1][0]
with open (sRNAmappings_file, 'r') as fh: 
  mappings = [x.rstrip('\n').split('\t') for x in fh.readlines()]

# newdata = []
# for m in mappings: 
  # srna = m[0]
  # freq = m[1]
  # strd = str(m[3])
  # if strd == '-': srna = getRevComp (srna)
  # if strd == '+': strd = 'p'
  # else: strd = 'n'
  # posch= m[5]
  # line = str(posch).zfill(4) + '.'*int(posch) + srna + '_' + freq + '_' + strd + '_' + str(len(srna))
  # newdata.append(line)
# newdata = sorted(newdata)
# for line in newdata:
  # print(line)
# print('....' + genseq)

#print(len(genseq))
d_pos_freq = {}
for i in range (1, len (genseq) + 1 ):
  d_pos_freq[i] = 0

for m in mappings:
  srna = m[0]
  freq = int(m[1])
  strd = m[3]
  posch= int(m[5])
  if strd == '-': posch = posch - 2
  if len(srna) == 21:
    d_pos_freq[posch] += freq
  
#for k, v in d_pos_freq.items():
#  print(k, v)


for i in range (2, len (genseq) + 1 - 21*7):
  n = 0
  k1 = d_pos_freq[i]
  if k1 > 0: n += 1
  k2 = d_pos_freq[i+21*1]
  if k2 > 0: n += 1
  k3 = d_pos_freq[i+21*2]
  if k3 > 0: n += 1
  k4 = d_pos_freq[i+21*3]
  if k4 > 0: n += 1
  k5 = d_pos_freq[i+21*4]
  if k5 > 0: n += 1
  k6 = d_pos_freq[i+21*5]
  if k6 > 0: n += 1
  k7 = d_pos_freq[i+21*6]
  if k7 > 0: n += 1
  k8 = d_pos_freq[i+21*7]
  if k8 > 0: n += 1
  sum1_8 = sum([k1, k2, k3, k4, k5, k6, k7, k8])
  print(i, sum1_8, n)


