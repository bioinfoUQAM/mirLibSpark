'''
Chao-Jung Wu
2018-09-29

positions are zero-based
'''

from __future__ import print_function
import math

def getRevComp (seq):
  intab = "ACGT"
  outab = "TGCA"
  #= v1 ================================
  #from string import maketrans
  #trantab = maketrans(intab, outab)
  #= v2 ================================
  trantab = str.maketrans(intab, outab)
  #=====================================
  n_seq = seq.translate(trantab)
  return n_seq[::-1]
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def calculate_phase_score (sum_of_k, n):
  #= y = pow (x, a) => x to the power of a
  #= y = math.log (x) => If base not defined, base is e as natural log. Take log of x with base e
  if n < 3: p = 0
  elif sum_of_k < 1: p = 0
  else:
    p = math.log(pow((1 + sum_of_k), n-2)) 
  return round(p, 3)

def display_alignment ():
  newdata = []
  for m in mappings: 
    srna = m[0]
    freq = m[1]
    strd = str(m[3])
    if strd == '-': srna = getRevComp (srna)
    if strd == '+': strd = 'p'
    else: strd = 'n'
    posch= m[5]
    line = str(posch).zfill(4) + '.'*int(posch) + srna + '_' + freq + '_' + strd + '_' + str(len(srna))
    newdata.append(line)
  newdata = sorted(newdata)
  for line in newdata:
    print(line)
  print('....' + genseq)

def get_d_Kbp_pos_freq (k):
  ''' 
  k: the length of wanted small RNAs 
  k = int or 'All', where All for len(smallRNA) > 17 bp
  '''
  d_Kbp_pos_freq = {}
  for i in range (1, len (genseq) + 1 ):
    d_Kbp_pos_freq[i] = 0

  for m in mappings:
    srna = m[0]
    freq = int(m[1])
    strd = m[3]
    posch= int(m[5])
    if strd == '-': posch = posch - 2

    if k == 'All':
      if len(srna) > 17:
        d_Kbp_pos_freq[posch] += freq
    elif len(srna) == k:
      d_Kbp_pos_freq[posch] += freq
      
  return d_Kbp_pos_freq 


def report_phase_scores (): 
  shift_position = 74 #= shoft to the middle of 4th cycle
  print( '\t'.join('position Reads_Sum_of_k n pos_mid4thcycle phaseScore'.split(' ')))
  for i in range (2, len (genseq) + 1 - 21*7):
    n = 0
    k1 = d_21bp_pos_freq[i]
    if k1 > 0: n += 1
    k2 = d_21bp_pos_freq[i+21*1]
    if k2 > 0: n += 1
    k3 = d_21bp_pos_freq[i+21*2]
    if k3 > 0: n += 1
    k4 = d_21bp_pos_freq[i+21*3]
    if k4 > 0: n += 1
    k5 = d_21bp_pos_freq[i+21*4]
    if k5 > 0: n += 1
    k6 = d_21bp_pos_freq[i+21*5]
    if k6 > 0: n += 1
    k7 = d_21bp_pos_freq[i+21*6]
    if k7 > 0: n += 1
    k8 = d_21bp_pos_freq[i+21*7]
    if k8 > 0: n += 1
    sum_of_k = sum([k1, k2, k3, k4, k5, k6, k7, k8])
    p = calculate_phase_score (sum_of_k, n)
    print( '\t'.join( str(j) for j in [i, sum_of_k, n, (shift_position + i), p])  )

####################################################################################################

genome_seq_file = '../dbs/JEAN/Genome/XLOC_126981_6D_REVCOMP.fasta'
sRNAmappings_file = 'demo_loci.txt' # bowtie result: seq	freq	nbloci	strand	chromo	poschr
with open (genome_seq_file, 'r') as fh: 
  genseq = [x.rstrip('\n').split('\t') for x in fh.readlines()][1][0]
with open (sRNAmappings_file, 'r') as fh: 
  mappings = [x.rstrip('\n').split('\t') for x in fh.readlines()]


def totalReads_ofeachLengthSmallRNA ():
  d_len_freq = {}
  for m in mappings:
    srna = m[0]
    lensrna = len(srna)
    freq = int(m[1])
    strd = m[3]
    posch= int(m[5])
    if strd == '-': posch = posch - 2
    if lensrna < 18: continue
    if lensrna not in d_len_freq.keys(): d_len_freq[lensrna] = freq
    else: d_len_freq[lensrna] += freq
  print('lenSmallRNA\ttotalFreq')
  for k in sorted (d_len_freq.keys()):
    v = d_len_freq[k]
    print(str(k) + '\t' + str(v))


def argonaute_nucleotide_preference_analysis ():
  d_21_5prime = {'A':0, 'T':0, 'C':0, 'G':0}
  d_21_3prime = {'A':0, 'T':0, 'C':0, 'G':0}
  d_22_5prime = {'A':0, 'T':0, 'C':0, 'G':0}
  d_22_3prime = {'A':0, 'T':0, 'C':0, 'G':0}
  for m in mappings:
    srna = m[0]
    lensrna = len(srna)
    freq = int(m[1])
    strd = m[3]
    posch= int(m[5])
    if lensrna == 21:
      first = srna[0]
      d_21_5prime[first] += freq
      last = srna[-1]
      d_21_3prime[last] += freq
    if lensrna == 22:
      first = srna[0]
      d_22_5prime[first] += freq
      last = srna[-1]
      d_22_3prime[last] += freq
  print('21bp_5prime\tTotalFreq')
  for k, v in sorted(d_21_5prime.items()): print(str(k) + '\t' + str(v))
  print('\n21bp_3prime\tTotalFreq')
  for k, v in sorted(d_21_3prime.items()): print(str(k) + '\t' + str(v))
  print('\n22bp_5prime\tTotalFreq')
  for k, v in sorted(d_22_5prime.items()): print(str(k) + '\t' + str(v))
  print('\n22bp_3prime\tTotalFreq')
  for k, v in sorted(d_22_3prime.items()): print(str(k) + '\t' + str(v))





#=================================
#=
#=  MAIN
#=
#=================================
#display_alignment ()

#= figure_1A_upperPannel
#d_18ormorebp_pos_freq = get_d_Kbp_pos_freq ('All')
#for k, v in d_18ormorebp_pos_freq.items(): print(str(k) + '\t' + str(v))

#= figure_1A_lowerPannel
#d_21bp_pos_freq = get_d_Kbp_pos_freq (21)
#report_phase_scores ()

#totalReads_ofeachLengthSmallRNA ()

argonaute_nucleotide_preference_analysis ()
