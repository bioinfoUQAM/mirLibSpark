'''
program: sequential.py
author: Chao-Jung Wu
date: 2017-10-11
version: 0.00.01
  
'''

from __future__ import print_function
import sys
import os.path
import time
#from os import listdir
#
import utils as ut
import mirLibRules as mru

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
      freq = data[1]
      dict_mirna[seq] = [freq]
  return dict_mirna

def filterFreq (limit_srna_freq, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna.items():
    if int(v[0]) < limit_srna_freq:
      del dict_mirna2[k]
  return dict_mirna2

def filterShort (limit_len, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k in dict_mirna.keys():
    if len(k) < limit_len:
      del dict_mirna2[k]
  return dict_mirna2

def filterdust (dict_mirna):
  #= echo $'>seq1\nCGTGGCTATGATAGCGATATTCGTTTTTTT' | dustmasker
  dict_mirna2 = dict_mirna.copy()
  for k in dict_mirna.keys():
      cmd = 'echo $\'>seq1\n' + k + '\' | dustmasker > dustmasker.tmp2'
      os.system(cmd)
      with open ('dustmasker.tmp2', 'r') as fh:
        data = fh.readlines()
        if len(data) == 2: 
          del dict_mirna2[k]
  return dict_mirna2
  
def bowtiemap (dict_mirna):
  bowtie_index = '/home/cloudera/workspace/mirLibHadoop/dbs/bowtie_index/a_thaliana_t10'
  for k in dict_mirna.keys():
    cmd = 'bowtie --mm -a -v 0 --suppress 1,5,6,7,8 -c ' + bowtie_index + ' '+ k + ' > bowtiemap.tmp2 2>/dev/null'
    os.system(cmd)
    locs = []
    with open ('bowtiemap.tmp2', 'r') as fh:
      for line in fh:
        data = line.rstrip('\n').split('\t') #= +	Chr2	1040947
        locs.append(data)
    dict_mirna[k].append(len(locs))
    dict_mirna[k].append(locs) #= v =[freq, 2, [[strand, chromo, pos],[strand, chromo, pos]]]

def filterMirnaFreq (limit_mrna_freq, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna.items():
    if int(v[0]) < limit_mrna_freq:
      del dict_mirna2[k]
  return dict_mirna2

def filterNbLoc (limit_nbloc, dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna.items():
    if v[1] == 0 or v[1] > limit_nbloc:
      del dict_mirna2[k]
  return dict_mirna2

def filterKnowNon (dict_mirna):
  dict_mirna2 = dict_mirna.copy()
  for k, v in dict_mirna.items():
    if not kn_obj.knFilterByCoor ([k, v]):
      del dict_mirna2[k]
  return dict_mirna2

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

#def extract_from_genome (strand, chromo, start_srna, len_srna):
#  if chromo not in extr_obj.genome.keys(): return [['-', -2]]
#  contig = extr_obj.genome[chromo]
#  prims = extr_obj.extract_precursors (contig, strand, start_srna, len_srna)
#  return prims

def fold1 (dict_mirna):
  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      prims = v[2][i][3]
      for j in range(len(prims)):
        cmd = 'echo ' + prims[j][0] + '| RNAfold  > rnafold.tmp2'
        os.system(cmd)
        with open ('rnafold.tmp2', 'r') as fh:
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
        cmd = 'perl eval_mircheck.pl \"' + folding + '\" ' + str(miRNA_start) + ' ' + str(miRNA_stop) + ' def > testcheck_tmp.txt'
        os.system(cmd)
        with open ('testcheck_tmp.txt', 'r') as fh:
          for line in fh:
            check = line.rstrip('\n').split('\t')
            if 'prime' in check[0]:
              dict_mirna[k][2][i][3][j].append(check[0])
              dict_mirna[k][2][i][3][j].append(check[1])
              dict_mirna[k][2][i][3][j].append(check[2])
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
        cmd = 'echo ' + preSeq + '| RNAfold  > rnafold2.tmp2'
        os.system(cmd)
        with open ('rnafold2.tmp2', 'r') as fh:
          fh.readline()
          for line in fh:
            fold = line.split(' (')[0]
            dict_mirna[k][2][i][3][j].append(fold)

def mirdupcheck (dict_mirna):
  #= cmd = 'java -jar ../lib/miRdup_1.4/miRdup.jar -v mirdupcheck.tmp2 -c ../lib/miRdup_1.4//model/thaliana.model -r /usr/local/bin/'
  pred_file = 'mirdupcheck.tmp2.thaliana.model.miRdup.txt'

  for k, v in dict_mirna.items():
    for i in range(len(v[2])):
      pres = v[2][i][3]
      for j in range(len(pres)):
        if len(pres[j]) == 0: continue 
        preSeq = pres[j][6]
        miseq = k
        fold = pres[j][8]
        #print(preSeq, miseq, fold)
  
        with open ('mirdupcheck.tmp2', 'w') as fhtmp:
          print('seqx\t' + miseq + '\t' + preSeq + '\t' + fold, file=fhtmp)

        cmd = 'java -jar ../lib/miRdup_1.4/miRdup.jar -v mirdupcheck.tmp2 -c ../lib/miRdup_1.4//model/thaliana.model -r /usr/local/bin/ 1>/dev/null'
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
  list_profile = []
  for k, v in dict_mirna.items():
    i = [k, v]
    list_profile.append(i)
  return list_profile

def filterProfile (dict_mirna):
  for k, v in dict_mirna.items():
    elem = [k, v]
    elem = profile_obj.computeProfileFrq(elem, dict_bowtie_chromo_strand)
    #print(elem)

    #profile_rdd = pre_vld_rdd.map(lambda e: profile_obj.computeProfileFrq(e, dict_bowtie_chromo_strand))\
     #                 .filter(lambda e: e[1][0] / float(e[1][5]) > 0.2)#\
                      #.persist()####################

        

#infile = 'test.txt'
infile = '/home/cloudera/Desktop/mirLibHadoop/input_storage/100.txt'
dict_mirna = readRaw (infile)
dict_mirna = filterFreq (limit_srna_freq, dict_mirna)
dict_mirna = filterShort (limit_len, dict_mirna)
dict_mirna = filterdust (dict_mirna)
bowtiemap (dict_mirna)
list_profile = dict_mirna_for_profile (dict_mirna) #########################
dict_mirna = filterMirnaFreq (limit_mrna_freq, dict_mirna)
dict_mirna = filterNbLoc (limit_nbloc, dict_mirna)
extractPri (dict_mirna)
fold1 (dict_mirna)
check (dict_mirna)
filterOneLoop (dict_mirna)
extractPre_fromPri (dict_mirna)
fold2 (dict_mirna)
mirdupcheck (dict_mirna)
#dict_bowtie_chromo_strand = profile_obj.get_bowtie_strandchromo_dict (list_profile)###########
#filterProfile (dict_mirna)


for k, v in dict_mirna.items():
  print(k,v)
  break




