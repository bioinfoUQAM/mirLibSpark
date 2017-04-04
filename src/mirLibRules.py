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

#=========================================================
class prog_remove_low_freq ():
  
  def __init__(self, vals_sep, kv_sep):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep

  def lowfreq_filter_rule(self, kv_arg):
    limit_freq = 6 # set up upper limit of frequency
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    freq = value.split(self.values_sep)[1]
    if int(freq) < limit_freq:
        return False
    return True  ## false data will automatically excluded in the new RDD



class prog_remove_short_length ():
  
  def __init__(self, vals_sep, kv_sep):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep

  def shortlen_filter_rule(self, kv_arg):
    limit_len = 14 # set up lower limit of the length of RNA sequence
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    RNAseq = value.split(self.values_sep)[0]
    if len(RNAseq) < limit_len:
        return False
    return True  ## false data will automatically excluded in the new RDD



class prog_remove_invalid_nbLocations ():
  
  def __init__(self, vals_sep, kv_sep):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep

  def nbLocations_filter_rule(self, kv_arg):
    limit_nbLocations = 5 # set up lower limit of nbLocation

    ## FACT ## len(kv_arg) = 1, it is a list, can not split

    #nbLoc_by_bowtie = len(kv_arg) - 1
    #if nbLoc_by_bowtie == 2:
    #    return True
    #return False
    
    #keyvalue = kv_arg.split(self.keyval_sep) # keyvalue = kv_arg.split(self.keyval_sep) AttributeError: 'list' object has no attribute 'split'
    #key = keyvalue[0]
    #value = keyvalue[1]

    if len(value) > limit_nbLocations + 1:
        return False
    if len(value) == 2:
        bowtie_result_1 = value[1]
        if bowtie_result_1[0] == 'NA':
            return False
    return True  ## false data will automatically excluded in the new RDD
#============================================================


class prog_dustmasker ():
  
  def __init__(self, vals_sep, kv_sep):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep

  def dmask_filter_rule(self, kv_arg):
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    sRNAseq = value.split(self.values_sep)[0]
    line1 = ['echo', '>seq1\n' + sRNAseq]
    line2 = ['dustmasker']
    p1 = sbp.Popen(line1, stdout=sbp.PIPE)
    p2 = sbp.Popen(line2, stdin=p1.stdout, stdout=sbp.PIPE)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate()[0].rstrip('\n')
    nblines = len(output.split('\n'))
    if nblines == 1:
        return True
    return False  ## false data will automatically excluded in the new RDD


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
