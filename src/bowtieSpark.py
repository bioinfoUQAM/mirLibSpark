'''
spark-submit bowtieSpark.py 2>/dev/null

Date: 2018-05-24
Author: Chao-Jung Wu
'''

from __future__ import print_function
import sys
import os.path
import time
from os import listdir
#
import utils as ut
import mirLibRules as mru


###################
#= Convert a fasta file into a key value file
#= (k, v) = (seq\tacc_num)
def covert_fasta_to_KV(infile, outfile):
  fh = open (infile, 'r')
  DATA = fh.readlines()
  fh.close()
  dict_sRNA = {}
  for i in xrange(0, len(DATA), 2):
      key = DATA[i].rstrip('\n')[1:]
      value = DATA[i+1].rstrip('\n')
      dict_sRNA[key] = value
  fh_out = open (outfile, 'w')
  for k, v in dict_sRNA.items():
      print (v + '\t' + k, file=fh_out)
  fh_out.close()

def rearrange_rule(kv_arg, kv_sep):
  tab = kv_arg.split(kv_sep)
  return (tab[0], tab[1])
####################

#= infile must be in seq\tint format
#infile = '/home/cloudera/Desktop/serialTGACv1.txt' 
infile = '/home/cloudera/Desktop/mirLibHadoop/input/fake_a.txt'
#outfile = '/home/cloudera/Desktop/mirLibHadoop/tmp/fake_c.kv.txt'
#covert_fasta_to_KV(infile, outfile)
#infile = outfile





appMaster, appName, mstrMemory, execMemory, execCores = 'local[*]', 'mirLibHadoop', '10g', '10g', '24'
b_index_path = '/home/cloudera/Desktop/mirLibHadoop/dbs/bowtie_index/a_thaliana_t10'
#b_index_path = '/home/cloudera/vm_share/180523_xloc_blowtie_index/bowtie_index_xloc/xlocIndex'

partition = 4


sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execCores)
bowtie_obj = mru.prog_bowtie(b_index_path)
bowtie_cmd, bowtie_env = bowtie_obj.Bowtie_pipe_cmd()


distFile_rdd = sc.textFile("file:///" + infile, partition).persist()
print(distFile_rdd.collect())

collapse_rdd = distFile_rdd.map(lambda line: rearrange_rule(line, '\t')).persist()
print(collapse_rdd.collect())


bowtie_rdd = distFile_rdd.map(lambda e: " "+e)\
                         .pipe(bowtie_cmd, bowtie_env)\
                         .map(bowtie_obj.bowtie_rearrange_map)\
                         .groupByKey()\
                         .map(lambda e: (e[0], [len(list(e[1])), list(e[1])]))\
                         .persist()
print(bowtie_rdd.collect())


bowFrq_rdd = bowtie_rdd.join(collapse_rdd)\
                       .map(bowtie_obj.bowtie_freq_rearrange_rule)\
                       .persist()


print(bowFrq_rdd.collect())