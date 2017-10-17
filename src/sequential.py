'''
program: sequential.py
author: Chao-Jung Wu
date: 2017-10-11
version: 0.00.02

100.txt before profile repeats: (1) 116m45s, (2) 96mins
  
'''

from __future__ import print_function
import sys
import os.path
import time
#from os import listdir
#
import utils as ut
import mirLibRules as mru

tmp_rep = '/home/cloudera/Desktop/mirLibHadoop/tmp/'

known_non = '../dbs/TAIR10_ncRNA_CDS.gff'
d_ncRNA_CDS = ut.get_nonMirna_coors (known_non)
kn_obj = mru.prog_knownNonMiRNA(d_ncRNA_CDS)
profile_obj = mru.prog_dominant_profile()

limit_srna_freq = 10
limit_len = 18 
limit_mrna_freq = 100
limit_nbloc = 15

genome_path = '/home/cloudera/workspace/mirLibHadoop/dbs/ATH/Genome/'
pri_l_flank = 20 #500
pri_r_flank = 160 #200
pre_flank = 30
extr_obj = mru.extract_precurosrs (genome_path, pri_l_flank, pri_r_flank, pre_flank)


def readRaw (infile):
  ''' seq\tfreq'''
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
      cmd = 'echo $\'>seq1\n' + k + '\' | dustmasker > ' + tmp_rep + 'dustmasker.tmp2'
      os.system(cmd)
      with open (tmp_rep + 'dustmasker.tmp2', 'r') as fh:
        data = fh.readlines()
        if len(data) == 2: 
          del dict_mirna[k]
  
def bowtiemap (dict_mirna):
  bowtie_index = '/home/cloudera/workspace/mirLibHadoop/dbs/bowtie_index/a_thaliana_t10'
  for k in dict_mirna.keys():
    cmd = 'bowtie --mm -a -v 0 --suppress 1,5,6,7,8 -c ' + bowtie_index + ' '+ k + ' > ' + tmp_rep + 'bowtiemap.tmp2 2>/dev/null'
    os.system(cmd)
    locs = []
    with open (tmp_rep + 'bowtiemap.tmp2', 'r') as fh:
      for line in fh:
        data = line.rstrip('\n').split('\t') #= +	Chr2	1040947
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
        cmd = 'echo ' + prims[j][0] + '| RNAfold  > ' + tmp_rep + 'rnafold.tmp2'
        os.system(cmd)
        with open (tmp_rep + 'rnafold.tmp2', 'r') as fh:
          fh.readline()
          for line in fh:
            fold = line.split(' (')[0]
            dict_mirna[k][2][i][3][j].append(fold)

def check (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      prims = v[2][i][3]
      for j in range(len(prims)):
        miRNA_start = int(prims[j][1])
        miRNA_stop = miRNA_start + len(k) -1
        folding = prims[j][2]
        cmd = 'perl eval_mircheck.pl \"' + folding + '\" ' + str(miRNA_start) + ' ' + str(miRNA_stop) + ' def > ' + tmp_rep + 'testcheck_tmp.txt'
        os.system(cmd)
        with open (tmp_rep + 'testcheck_tmp.txt', 'r') as fh:
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
        preSeq = pres[j][6]
        cmd = 'echo ' + preSeq + '| RNAfold  > ' + tmp_rep + 'rnafold2.tmp2'
        os.system(cmd)
        with open (tmp_rep + 'rnafold2.tmp2', 'r') as fh:
          fh.readline()
          for line in fh:
            fold = line.split(' (')[0]
            dict_mirna[k][2][i][3][j].append(fold)

def mirdupcheck (dict_mirna):
  #= cmd = 'java -jar ../lib/miRdup_1.4/miRdup.jar -v mirdupcheck.tmp2 -c ../lib/miRdup_1.4//model/thaliana.model -r /usr/local/bin/'
  pred_file = tmp_rep + 'mirdupcheck.tmp2.thaliana.model.miRdup.txt'

  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      pres = v[2][i][3]
      for j in range(len(pres)):
        if len(pres[j]) == 0: continue 
        preSeq = pres[j][6]
        miseq = k
        fold = pres[j][8]
        #print(preSeq, miseq, fold)
  
        with open (tmp_rep + 'mirdupcheck.tmp2', 'w') as fhtmp:
          print('seqx\t' + miseq + '\t' + preSeq + '\t' + fold, file=fhtmp)

        cmd = 'java -jar ../lib/miRdup_1.4/miRdup.jar -v ' + tmp_rep + 'mirdupcheck.tmp2 -c ../lib/miRdup_1.4//model/thaliana.model -r /usr/local/bin/ 1>/dev/null'
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
  '''
['GCTCACTGCTCTTTCTGTCAGA', ['500', 1, [['-', 'Chr2', '10676469', [['TTCTCATCGTTTCTTGTTTTCTTTGTTTCATCTTGTAGATCTCTGAAGTTGGACTAATTGTGAATGAAAGAGTTGGGACAAGAGAAACGCAAAGAAACTGACAGAAGAGAGTGAGCACACAAAGGCAATTTGCATATCATTGCACTTGCTTCTCTTGCGTGCTCACTGCTCTTTCTGTCAGATTCCGGTGCTGATCTCTTTG', 160, '.((.(((((........((((((((((((.((((((....((((...(((..((.....)).)))...)))).....)))))))))..)))))))))(((((((((((((((((((((.(((((((((..((((......)))).))))))...))).))))))))).))).)))))))))....)))))..))........', '3prime', '96', '181', 'AAAGAGTTGGGACAAGAGAAACGCAAAGAAACTGACAGAAGAGAGTGAGCACACAAAGGCAATTTGCATATCATTGCACTTGCTTCTCTTGCGTGCTCACTGCTCTTTCTGTCAGATTCCGGTGCTGATCTCTTTG', 94, '.............((((((...(((..(((.(((((((((((((((((((((.(((((((((..((((......)))).))))))...))).))))))))).))).))))))))).)))...)))...))))))..', 'true', '0.72'], []]]]]]

elem = (seq, [frq, nbloc, [bowtie], [pri_miRNA], [pre_miRNA]])
 '''
  dict_mirna2 = dict_mirna.copy()

  for k, v in dict_mirna2.items():
    frq = int(v[0])
    nbLoc = v[1]
    bowtie = v[2][0][0:2]
    bowtie.append(int(v[2][0][2]))
    items = v[2][0][3]
    for i in items:
      if len(i) == 0: continue
      prim = i[0:6]
      prem = i[6:11]
      elem = [k, [int(frq), nbLoc, bowtie, prim, prem]]

      elem = profile_obj.computeProfileFrq(elem, dict_bowtie_chromo_strand)
      #print(elem)
      

    #profile_rdd = pre_vld_rdd.map(lambda e: profile_obj.computeProfileFrq(e, dict_bowtie_chromo_strand))\
     #                 .filter(lambda e: e[1][0] / float(e[1][5]) > 0.2)#\
                      #.persist()####################

        

infile = 'test.txt'
#infile = '/home/cloudera/Desktop/mirLibHadoop/input/100.txt'
dict_mirna = readRaw (infile)
filterFreq (limit_srna_freq, dict_mirna)
filterShort (limit_len, dict_mirna)
filterdust (dict_mirna)
bowtiemap (dict_mirna)
#d2 = dict_mirna.copy()
#bowtie_collect = dict_mirna_for_profile (d2) #########################
#dict_bowtie_chromo_strand = profile_obj.get_bowtie_strandchromo_dict (bowtie_collect)
#for k, v in dict_bowtie_chromo_strand.items():
#  print(k,v)
filterMirnaFreq (limit_mrna_freq, dict_mirna)
filterNbLoc (limit_nbloc, dict_mirna)
extractPri (dict_mirna)
fold1 (dict_mirna)
check (dict_mirna)
filterOneLoop (dict_mirna)
extractPre_fromPri (dict_mirna)
fold2 (dict_mirna)
mirdupcheck (dict_mirna)


#filterProfile (dict_mirna)


for k, v in dict_mirna.items():
  print(k,v)
  #break




