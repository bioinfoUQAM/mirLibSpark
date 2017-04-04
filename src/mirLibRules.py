'''
program: mirLibRules.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 0.00.01

Le programme 

'''

import os, sys
import subprocess as sbp


class prog_bowtie ():
  
  def __init__(self, vals_sep, kv_sep, b_index):
    self.bowtie_index = b_index
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep
    
    
  def run_bowtie(self, seq):
    append_values = []
    
    # self.cmd = 'bowtie --mm -a -v 0 --suppress 1,5,6,7,8 -c ' + self.bowtie_index + ' '+ seq  # shell=True
    self.cmd = ['bowtie', '--mm', '-a', '-v', '0', '--suppress', '1,5,6,7,8', '-c', self.bowtie_index, seq] # shell=False
    
    sproc = sbp.Popen(self.cmd, stdout=sbp.PIPE, shell=False)
    bsout = sproc.communicate()[0]
    bwout = bsout.decode("ascii").rstrip('\n')
    
    if bwout :
      bwList = bwout.split('\n')
      for line in bwList:
        append_value = line.rstrip('\n').split('\t')
        append_values += [append_value]
      
    else :
      append_values = [['NA', 'NA', 'NA']]
    
    return append_values
  
  def Bowtie_map_rule(self, kv_arg):
    
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    sRNAseq = value.split(self.values_sep)[0]
    
    append_value = self.run_bowtie(sRNAseq)
    
    kv2 = [key, [value] + append_value]
    return kv2

if __name__ == '__main__' :
   
   values_sep = "<>"
   keyval_sep = "::"
   b_index = "a_thaliana_t10"
   
   bowtie = prog_bowtie(values_sep, keyval_sep, b_index)
   print(bowtie.run_bowtie('CAAAGACTCATATGGACTTTGG'))

'''
http://stackoverflow.com/questions/30010939/python-subprocess-popen-error-no-such-file-or-directory
'''