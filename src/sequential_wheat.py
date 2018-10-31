'''

program: sequential.py
author: Chao-Jung Wu
date: 2017-10-11
version: 0.00.03


update: 2018-10-28
Prepare to process 11w2013_t2_1.fasta
  

'''



from __future__ import print_function
import sys
import os
import os.path
import time
import datetime



#
import utils as ut
import mirLibRules as mru
import uuid



uid = uuid.uuid4().hex
adapter='TGGAATTCTCGGGTGCCAAGGAACTC'
#adapter='none'

input_type='fasta'





rep_tmp = '../tmp/'


known_non = '../dbs/ATH_TAIR10/TAIR10_ncRNA_CDS.gff'
d_ncRNA_CDS = ut.get_nonMirna_coors (known_non)
kn_obj = mru.prog_knownNonMiRNA(d_ncRNA_CDS)
profile_obj = mru.prog_dominant_profile()


limit_srna_freq = 10
limit_len = 18 
limit_mrna_freq = 100
limit_nbloc = 15

genome_path = '../dbs/WHEAT_IWGSC/Genome/'
bowtie_index = '../dbs/WHEAT_IWGSC/bowtie_index/'

pri_l_flank = 500 #20
pri_r_flank = 200 #160
pre_flank = 10

genome = ut.getGenome(genome_path, ".fa")
extr_obj = mru.extract_precurosrs (genome, pri_l_flank, pri_r_flank, pre_flank)

RNAfold_path = ut.find_RNAfold_path ()


def find_str(s, char):
    ''' zero based indexing '''
    index = 0
    if char in s:
        c = char[0]
        for ch in s:
            if ch == c:
                if s[index:index+len(char)] == char:
                    return index
            index += 1
    return -1000

def trim_adapter (seq, ad):
  '''
  example:  adapter ad =                  TGGAATTCTCGGGTGCCAAGGAACTC
            seq =        NTACCGATCTGAGCCATTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACN
            return =     NTACCGATCTGAGCCAT
  updated: 2018-10
  '''
  while len(ad) > 6:
    len_ad = len(ad)
    pos = find_str(seq, ad)
    if pos > 0: return seq[:pos]
    else: ad = ad[:-1]
  return seq


def readFasta(infile):
  dict_mirna = {}
  with open (infile, 'r') as fh:
    for line in fh:
      if not line.startswith ('>'):
        seq = line.rstrip('\n')
        if not adapter == 'none': seq = trim_adapter (seq, adapter)
        if not seq in dict_mirna: dict_mirna[seq] = [1]
        else: dict_mirna[seq][0] += 1
  return dict_mirna


def readRaw (infile):
  ''' seq\tfreq '''
  dict_mirna = {}
  with open (infile, 'r') as fh:
    for line in fh:
      data = line.rstrip('\n').split('\t')
      seq = data[0]
      freq = int(data[1])
      dict_mirna[seq] = [freq]
  return dict_mirna



def filterFreq (limit_srna_freq, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna2.items():
    if int(v[0]) < limit_srna_freq:
      del dict_mirna[k]



def filterShort (limit_len, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k in dict_mirna2.keys():
    if len(k) < limit_len:
      del dict_mirna[k]



def filterdust (dict_mirna):

  #= echo $'>seq1\nCGTGGCTATGATAGCGATATTCGTTTTTTT' | dustmasker
  dict_mirna2 = dict_mirna.copy()
  for k in dict_mirna2.keys():
      tmpfile = rep_tmp + uid + '_dustmasker.tmp2'
      cmd = 'echo $\'>seq1\n' + k + '\' | dustmasker > ' + tmpfile
      os.system(cmd)
      with open (tmpfile, 'r') as fh:
        data = fh.readlines()
        if len(data) == 2: 
          del dict_mirna[k]
  
def bowtiemap (dict_mirna):

  for k in dict_mirna.keys():
    tmpfile = rep_tmp + uid + '_bowtiemap.tmp2'
    cmd = 'bowtie --mm -a -v 0 --suppress 1,5,6,7,8 -c ' + bowtie_index + ' ' + k + ' >' + tmpfile + ' 2>/dev/null'
    os.system(cmd)
    locs = []
    with open (tmpfile, 'r') as fh:
      for line in fh:
        data = line.rstrip('\n').split('\t') #= +	Chr2	1040947
        strand = data[0]
        chromo = data[1]
        posSeq = int(data[2])
        data = [strand, chromo, posSeq]
        locs.append(data)
    dict_mirna[k].append(len(locs))
    dict_mirna[k].append(locs) #= v =[freq, 2, [[strand, chromo, pos],[strand, chromo, pos]]]

def filterMirnaFreq (limit_mrna_freq, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna2.items():
    if int(v[0]) < limit_mrna_freq:
      del dict_mirna[k]

def filterNbLoc (limit_nbloc, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna2.items():
    if v[1] == 0 or v[1] > limit_nbloc:
      del dict_mirna[k]

def filterKnowNon (dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna.items():
    if not kn_obj.knFilterByCoor ([k, v]):
      del dict_mirna2[k]

def extractPri (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      bow = v[2][i]
      strand = bow[0]
      chromo = bow[1]
      start_srna = int(bow[2])
      len_srna = len(k)
      if chromo not in extr_obj.genome.keys(): 
        prims = [['-', -2]]
      else:
        contig = extr_obj.genome[chromo]
        prims = extr_obj.extract_precursors (contig, strand, start_srna, len_srna)
      dict_mirna[k][2][i].append(prims)

def fold1 (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      prims = v[2][i][3]
      for j in range(len(prims)):
        tmpfile = rep_tmp + uid + '_rnafold.tmp2'
        cmd = 'echo ' + prims[j][0] + '| RNAfold  > ' + tmpfile
        os.system(cmd)
        with open (tmpfile, 'r') as fh:
          fh.readline()
          for line in fh:
            fold = line.split(' (')[0]
            dict_mirna[k][2][i][3][j].append(fold)

def check (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      prims = v[2][i][3]
      for j in range(len(prims)):
        tmpfile = rep_tmp + uid + '_testcheck_tmp.txt'
        miRNA_start = int(prims[j][1])
        miRNA_stop = miRNA_start + len(k) -1
        folding = prims[j][2]
        cmd = 'perl eval_mircheck.pl \"' + folding + '\" ' + str(miRNA_start) + ' ' + str(miRNA_stop) + ' def > ' + tmpfile
        os.system(cmd)
        with open (tmpfile, 'r') as fh:
          for line in fh:
            check = line.rstrip('\n').split('\t')
            if 'prime' in check[0]:
              dict_mirna[k][2][i][3][j].append(check[0])
              dict_mirna[k][2][i][3][j].append(int(check[1]))
              dict_mirna[k][2][i][3][j].append(int(check[2]))
            else:
              del dict_mirna[k][2][i][3][j][:]

def filterOneLoop (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      prims = v[2][i][3]
      for j in range(len(prims)):
        if len(prims[j]) == 0: continue
        folding = prims[j][2]
        s = int(prims[j][4])
        e = int(prims[j][5])+1
        if not ut.containsOnlyOneLoop (folding[s:e]):
          del dict_mirna[k][2][i][3][j][:]

def extractPre_fromPri (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      prims = v[2][i][3]
      for j in range(len(prims)):
        if len(prims[j]) == 0: continue
        prim = prims[j]
        priSeq = prim[0]
        posMir = int(prim[1])
        fback_start = int(prim[4])
        fback_stop = int(prim[5])
        data = extr_obj.extract_sub_seq (priSeq, posMir, fback_start, fback_stop)
        dict_mirna[k][2][i][3][j].append(data[0])
        dict_mirna[k][2][i][3][j].append(data[1])

def fold2 (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      pres = v[2][i][3]
      for j in range(len(pres)):
        if len(pres[j]) == 0: continue
        tmpfile = rep_tmp + uid + '_rnafold2.tmp2'
        preSeq = pres[j][6]
        cmd = 'echo ' + preSeq + '| RNAfold  > ' + tmpfile
        os.system(cmd)
        with open (tmpfile, 'r') as fh:
          fh.readline()
          for line in fh:
            fold = line.split(' (')[0]
            dict_mirna[k][2][i][3][j].append(fold)

def mirdupcheck (dict_mirna):
  #= cmd = 'java -jar ../lib/miRdup_1.4/miRdup.jar -v mirdupcheck.tmp2 -c ../lib/miRdup_1.4//model/thaliana.model -r /usr/local/bin/'
  pred_file = rep_tmp + 'mirdupcheck.tmp2.thaliana.model.miRdup.txt'

  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      pres = v[2][i][3]
      for j in range(len(pres)):
        if len(pres[j]) == 0: continue 
        preSeq = pres[j][6]
        miseq = k
        fold = pres[j][8]
        #print(preSeq, miseq, fold)
        tmpfile = rep_tmp + uid + '_mirdupcheck.tmp2'
  
        with open (tmpfile, 'w') as fhtmp:
          print('seqx\t' + miseq + '\t' + preSeq + '\t' + fold, file=fhtmp)

        cmd = 'java -jar ../lib/miRdup_1.4/miRdup.jar -v ' + tmpfile + ' -c ../lib/miRdup_1.4//model/thaliana.model -r ' + RNAfold_path + ' 1>/dev/null'
        os.system(cmd)

        with open (pred_file, 'r') as fh_tmp :
          for line in fh_tmp :
            if line.startswith("#PR") :
              mirdup_pred = line.rstrip("\n").split("\t")[2]
            elif line.startswith("#SC") :
              mirdup_score = '%.2f' % round(float(line.rstrip("\n").split("\t")[2]), 2)
        dict_mirna[k][2][i][3][j].append(mirdup_pred)
        dict_mirna[k][2][i][3][j].append(mirdup_score)

def dict_mirna_for_profile (dict_mirna):
  ''' 'elem : (seq, [frq, nbloc, [bowties]])  '''
  list_profile = []
  for k, v in dict_mirna.items():
    i = [k, v]
    list_profile.append(i)
  return list_profile

def filterProfile (dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna2.items():
    frq = v[0]
    nbLoc = v[1]
    locs = v[2]
    for l in range(len(locs)):
      loc = locs[l]
      bowtie = loc[0:3]
      items = loc[3]
      for i in range(len(items)):
        item = items[i]
        if len(item) == 0: continue
        prim = item[0:6]
        prem = item[6:11]
        elem = [k, [frq, nbLoc, bowtie, prim, prem]]
        elem = profile_obj.computeProfileFrq(elem, dict_bowtie_chromo_strand)
        totalFrq = elem[1][5]
        ratio = frq / (float (totalFrq) + 0.1)
        if ratio > 0.2:
          dict_mirna[k][2][l][3][i].append(ratio)
        else:
          del dict_mirna[k][2][l][3][i][:]
  

def get_bowtie_strandchromo_dict (bowtie_rdd_collect):
    '''elem : (seq, [frq, nbloc, [bowties]])
    '''
    dict_bowtie_chromo_strand = {}
    
    for elem in bowtie_rdd_collect :
      bowties = elem[1][2]
      
      for bowtie in bowties :
        #= concatenate chromosome (bowtie[1]) and strand (bowtie[0])
        chromo_strand = bowtie[1] + bowtie[0]
        
        if chromo_strand not in dict_bowtie_chromo_strand.keys():
          dict_bowtie_chromo_strand[chromo_strand] = []
        
        dict_bowtie_chromo_strand[chromo_strand].append(elem)
    
    return dict_bowtie_chromo_strand

def createBowFrqDict (dict_mirna):
  ''' note that it must use copy module to make a copy for this operation, otherwise, the return dict_bowtie_chromo_strand will be mutabel with dict_mirna 
      By using the copied d2, it dissociates the mutation link. This might be a particular inconvenience in python 2.7
  '''
  import copy
  d2 = copy.deepcopy(dict_mirna)
  bowtie_collect = dict_mirna_for_profile (d2) 
  dict_bowtie_chromo_strand = get_bowtie_strandchromo_dict (bowtie_collect)   
  return dict_bowtie_chromo_strand

def keepTrue(dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  count = 0 #= count distinct chromo strand posChr valide precursor
  for k, v in dict_mirna2.items():
    TRUE = 0
    frq = v[0]
    nbLoc = v[1]
    locs = v[2]
    for l in range(len(locs)):
      loc = locs[l]
      bowtie = loc[0:3]
      items = loc[3]
      for i in range(len(items)):
        item = items[i]
        if not len(item) == 0: 
          TRUE = 1
      if TRUE == 1: count += 1
    if TRUE == 0:
      del dict_mirna[k]
  return count





infile = '../input_samples/11w2013_t2_1.fasta'



time_a = datetime.datetime.now()


#dict_mirna = readRaw (infile)

print(time_a, 'start')
dict_mirna = readFasta (infile)
for k, v in dict_mirna.items(): 
    print('1st item after readFasta (infile): ', k, v)
    break
print(datetime.datetime.now(), 'dict_mirna = readFasta (infile)')
filterFreq (limit_srna_freq, dict_mirna)
print(datetime.datetime.now(), 'filterFreq (limit_srna_freq, dict_mirna)')
filterShort (limit_len, dict_mirna)
print(datetime.datetime.now(), 'filterShort (limit_len, dict_mirna)')
filterdust (dict_mirna)
print(datetime.datetime.now(), 'filterdust (dict_mirna)')
bowtiemap (dict_mirna)
print(datetime.datetime.now(), 'bowtiemap (dict_mirna)')
dict_bowtie_chromo_strand = createBowFrqDict (dict_mirna)
print(datetime.datetime.now(), 'dict_bowtie_chromo_strand = createBowFrqDict (dict_mirna)')
filterMirnaFreq (limit_mrna_freq, dict_mirna)
print(datetime.datetime.now(), 'filterMirnaFreq (limit_mrna_freq, dict_mirna)')
filterNbLoc (limit_nbloc, dict_mirna)
print(datetime.datetime.now(), 'filterNbLoc (limit_nbloc, dict_mirna)')
extractPri (dict_mirna)
print(datetime.datetime.now(), 'extractPri (dict_mirna)')
fold1 (dict_mirna)
print(datetime.datetime.now(), 'fold1 (dict_mirna)')
check (dict_mirna)
print(datetime.datetime.now(), 'check (dict_mirna)')
filterOneLoop (dict_mirna)
print(datetime.datetime.now(), 'filterOneLoop (dict_mirna)')
extractPre_fromPri (dict_mirna)
print(datetime.datetime.now(), 'extractPre_fromPri (dict_mirna)')
fold2 (dict_mirna)
print(datetime.datetime.now(), 'fold2 (dict_mirna)')
mirdupcheck (dict_mirna)
print(datetime.datetime.now(), 'mirdupcheck (dict_mirna)')
filterProfile (dict_mirna)
print(datetime.datetime.now(), 'filterProfile (dict_mirna)')
nbDistinctTrue = keepTrue(dict_mirna)
print(datetime.datetime.now(), 'nbDistinctTrue = keepTrue(dict_mirna)')


for k, v in dict_mirna.items(): print(k, v)

print('nbDistinctTrue: ', nbDistinctTrue)
print('nbUniqueMiRNA: ', len(dict_mirna))

time_b = datetime.datetime.now()
print(time_b, 'finish time')
print('total running time: ', time_b - time_a)


'''
#= final report one sequence example
AAGACTCATGAATTCATACGACAT [238, 1, [['+', '3', 15682287, [[], ['TGTTTTGCTTATGTTATGCTTATGTTTTTCTTATTTTTGTTCGATTATTTTATTTTTCTTAGTTTATTAGTCAGGACACATCAAAAAATTGCTACCATATTGTGACAATCTAATGAACACATCAGAAAGCCAAAGAAACTGGAGAAATAAGTGTGTGTTGGGCATAAAAAAAGAATAGAGATTTGGTGTTGATATGGTTTAAGACTCATGAATTCATACGACATGACCCTCCATAAAATATGTAACTACCCCAATCCTTTCTAAACATGACTCTGGTGTTGTAGAACTTATGATCGTCTACTTTTACTACTACGTCTTGAACGGGCTTGATTATGTTGTTATATACGATCTCCAAGTTGCGAAAGTAGCTAACTTTGTCAAAGCCTTCATCGGAGAAGTGACCAGATCCCATCTGTGTAATTGTGTGTCAGCCAAAACTTTGCGAATTCACAACCTCACCACCCCATTGCACGTCGTATATATAAGCTAAATTTGGAAATTTGCAGTTAACCAATATCCAACCAATGAATGGTTGGATCCCAAACTTAGCCACCAGTTCCCAGACTCTTGATCCTATATAATTGTTTATCAACGATGATATTTAACTTTTAAATTAGAAATGTAAAATTCATTAACTAAAACAAAAGATACAATCTCAAATATTTTTCCTCGGATGTAGAGACGGAGTAATGGTGCCTCCAACAACAATGTTGTTACTTGTTTG', 200, '.......(((...(((((((((.......((((((((..((((.((.((((..((((((..(((((((((....(.((((....................)))).)...))))))))).....))))))..)))).)).)))))))))))).......)))))))))...)))....(((((.(((((.(((...(((...((((((((((.(((.(((((((.......(((.....(((((.....................)))))....))))))))))))).))))))..))))))).))).))))).)))))(((((((..((((.(((((((((((((((((.....((((((....))))))......((((....((((......)))))))).((((...))))(((((((((((.((..((.........))..)).)))).............)))))))))))))))))...(((((.(((((....(((........))).(((((((((.....))))))))).))))).)))))(((((...(((.(.((((..((((.((((.((((((....)))))).))))....(((((...((((.((((.......))))..))))...)))))........................))))..)))).))).)...))))).........))))))).))))))))))).', '5prime', 200, 299, 'AAGAATAGAGATTTGGTGTTGATATGGTTTAAGACTCATGAATTCATACGACATGACCCTCCATAAAATATGTAACTACCCCAATCCTTTCTAAACATGACTCTGGTGTTGTAGAACTTATGATCGTCTACTTTTACTACTACGTCTTGAACGGGCTTGA', 30, '.......(((((.(((((.(((...(((...((((((((((.(((.(((((((.......(((.....(((((.....................)))))....))))))))))))).))))))..))))))).))).))))).)))))............', 'true', '0.67', 2380.0]]]]]


seq [freq, nbLoc,[[strand, chromo, posChr, [[], [ priSeq, 200, priFold, '5prime', 200, 299, preSeq, 30, preFold, 'true', '0.67', 2380.0]]]]]

'''



