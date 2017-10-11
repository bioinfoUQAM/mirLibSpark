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
#import utils as ut
#import mirLibRules as mru

limit_srna_freq = 10
limit_len = 18 

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

def dmask (dict_mirna):
  #= echo '>seqMir\n' + sRNAseq + '| dustmasker'
  

infile = 'test.txt'
dict_mirna = readRaw (infile)
dict_mirna = filterFreq (limit_srna_freq, dict_mirna)
dict_mirna = filterShort(limit_len, dict_mirna)

for k, v in dict_mirna.items():
  print(k,v)


