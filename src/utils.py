'''
program: utils.py
author: Chao-Jung Wu
author: M.A.Remita
date: 2017-03-25
update: 2019-10-16
version: 1.00.02

support python3 print syntax


#= 2018-10-28 note:
find a time to refactor diff analysis and enrichment analysis, make them into objects
'''
from __future__ import print_function
import os
import os.path
#import re
#import subprocess
#import sys
#from os import listdir

pythonV = 2


def transpose_txt(infile, outfile):
    with open(infile, 'r') as f:
        lis = [x.rstrip('\n').split('\t') for x in f]
    fho = open (outfile, 'w')
    for x in zip(*lis):
        for y in x:
            print(y+'\t', end='', file=fho)
        print('', file=fho)

def makedirs_reps (reps):
  for rep in reps:
    if not os.path.exists(rep):
      os.makedirs(rep)

# source : https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python/23728630#23728630
def randomStrGen (n):
  import string, random
  return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(n))

def roundup(x): 
  return x if x % 10 == 0 else x + 10 - x % 10

# Configure a spark context
def pyspark_configuration(appName, masterMemory, heartbeat):
  from pyspark import SparkConf, SparkContext
  myConf = SparkConf()
  myConf.setAppName(appName)  #= 'mirLibSpark'
  myConf.set("spark.driver.maxResultSize", '1500M') #= default 1g
  myConf.set("spark.driver.memory", masterMemory)
  timeout = heartbeat * 12
  myConf.set('spark.network.timeout', str(timeout) + 's')
  myConf.set('spark.executor.heartbeatInterval', str(heartbeat) + 's')
  return SparkContext(conf = myConf)

def convertTOhadoop(rfile, hdfsFile):
  '''
  Convert a file to hadoop file
  defunct
  '''
  print('if pre-existing in hdfs, the file would be deleted before the re-distribution of a new file with the same name.\n')
  os.system('hadoop fs -rm ' + hdfsFile) # force delete any pre-existing file in hdfs with the same name.
  os.system('hadoop fs -copyFromLocal ' + rfile + ' ' + hdfsFile)

def convert_fastq_file_to_KeyValue(infile, rep_tmp, inBasename):
  '''
  1: @SEQ_ID
  2: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  3: +
  4: !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
  outfile:
  seq\tquality(4th line)
  '''
  fh = open (infile, 'r')
  outfile = rep_tmp + inBasename + '.fastqKV'
  fh_out = open (outfile, 'w')
  i = 1
  for line in fh:
    if i == 1 or i == 3: 
      i += 1
      continue
    elif i == 2:
      seq = line.rstrip('\n')
      i += 1
    else:
      quality = line.rstrip('\n')
      print (seq + '\t' + quality, file=fh_out)
      i = 1
      continue
  fh.close()
  fh_out.close()
  return outfile

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


#===========================================================
#===========================================================
#===========================================================
#===========================================================
def stringTrans (v, intab, outab):
  if v == 2:
    from string import maketrans
    trantab = maketrans(intab, outab)
  if v == 3:
    trantab = str.maketrans(intab, outab)
  return trantab

def getRevComp (seq):
  intab, outab = "ACGT", "TGCA"
  trantab = stringTrans (pythonV, intab, outab)
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def tr_T_U (seq):
  intab, outab = "T", "U"
  trantab = stringTrans (pythonV, intab, outab)
  return seq.translate(trantab)

def tr_U_T (seq):
  intab, outab = "U", "T"
  trantab = stringTrans (pythonV, intab, outab)
  return seq.translate(trantab)
#===========================================================
#===========================================================
#===========================================================

def getChromosomeName (file):
  desc = ""
  with open(file, "r") as fh:
    for line in fh :
      if line.startswith(">"):
        desc = line
        break
  fh.close()
  return desc.split()[0][1:]

def getFastaSeq (file):
  seq = ""
  with open(file, "r") as fh:
    for line in fh :
      if not line.startswith(">"):
        seq = seq + line.rstrip("\n")
  fh.close()
  return seq

def getGenome (genome_path, file_ext, chromosomeName='All'):
  '''
  take in either All as the entire genome, or one chromosome at a time.
  '''
  genome = dict()
  if chromosomeName == 'All':
    files = [each for each in os.listdir(genome_path) if each.endswith(file_ext)]
  else: files = [ chromosomeName + file_ext ]

  for namefile in files :
    file = genome_path + namefile
    chr = getChromosomeName(file)
    sequence = getFastaSeq(file)
    genome[chr] = sequence
  return genome

def readParam (paramfile, sep = '='):
  paramDict = {}
  fh = open (paramfile, 'r')
  DATA = fh.readlines()
  fh.close()
  for line in DATA:
    if line.startswith('message'): 
      msg = line.rstrip('\r\n')[8:].split('\\n')
      paramDict['message'] = '\n'.join(msg)
    elif not line.startswith("#"):
      data = line.rstrip('\r\n').split(sep)
      paramDict[data[0]] = data[1]
  return paramDict

def writeTimeLibToFile (timeDict, outfile, appId, paramDict):
  import datetime
  
  totalTimeSec = reduce((lambda x,y : x + y), timeDict.values())
  totalTimeHMS = str(datetime.timedelta(seconds=totalTimeSec))
  totalTimeSec = format(totalTimeSec, '.3f')
  
  fh_out = open (outfile, 'w')
  print(datetime.datetime.now(), file=fh_out)
  print("# Application ID " + appId +"\n", file=fh_out)
  print("Total \t"+totalTimeHMS+ "\t" + totalTimeSec, file=fh_out)

  for lib in timeDict :
    timeLibSec = timeDict[lib]
    timeLibHMS = str(datetime.timedelta(seconds=timeLibSec))
    timeLibSec = format(timeLibSec, '.3f')
    print("Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec, file=fh_out)
    print ("Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec) #= stdout
  
  print("\n# SPARK configuration:", file=fh_out)

  for key in paramDict :
    if key.startswith("sc_"):
      print("# " + key + ": " + paramDict[key], file=fh_out)

  print("\n# MirLibSpark configuration:", file=fh_out)

  for key in sorted(paramDict.keys()) :
    if not key.startswith("sc_"):
      print("# " + key + ": " + paramDict[key], file=fh_out)
 
  fh_out.close()
