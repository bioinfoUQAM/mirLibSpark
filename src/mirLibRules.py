'''
program: mirLibRules.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 0.00.01

Le programme 

'''

import subprocess as sbp


class prog_remove_low_freq ():
  
  def __init__(self, vals_sep, kv_sep, limit_freq):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep
    self.limit_freq = limit_freq

  def lowfreq_filter_rule(self, kv_arg):
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    freq = value.split(self.values_sep)[1]
    if int(freq) < self.limit_freq:
        return False
    return True


class prog_remove_short_length ():
  
  def __init__(self, vals_sep, kv_sep, limit_len):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep
    self.limit_len = limit_len

  def shortlen_filter_rule(self, kv_arg):
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    RNAseq = value.split(self.values_sep)[0]
    if len(RNAseq) < self.limit_len:
        return False
    return True


class prog_remove_invalid_nbLocations ():
  
  def __init__(self, vals_sep, kv_sep, limit_nbLoc):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep
    self.limit_nbLoc = limit_nbLoc

  def nbLocations_filter_rule(self, kv_arg):
    # kv_arg = [u'000000001', [u'ATACGATCAACTAGAATGACAATT<>20', [u'-', u'Chr4', u'11833108']]]
    key = kv_arg[0] # 000000001
    values = kv_arg[1] # [u'ATACGATCAACTAGAATGACAATT<>20', [u'-', u'Chr4', u'11833108']]
    nbLoc = len(values) - 1
    if nbLoc > self.limit_nbLoc:
      return False
    bowtie_location_result_1 = values[1]
    if bowtie_location_result_1 == ['NA', 'NA', 'NA']: # remove zero location mapped by this sequence
      return False
    return True


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
    return False  ## false data will be automatically excluded in the new RDD


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

'''
#=== partition ========================================================================
def portable_hash(x):
    """
    This function returns consistent hash code for builtin types, especially
    for None and tuple with None.

    The algorithm is similar to that one used by CPython 2.7

    >>> portable_hash(None)
    0
    >>> portable_hash((None, 1)) & 0xffffffff
    219750521
    """
    if sys.version >= '3.3' and 'PYTHONHASHSEED' not in os.environ:
        raise Exception("Randomness of hash of string should be disabled via PYTHONHASHSEED")

    if x is None:
        return 0
    if isinstance(x, tuple):
        h = 0x345678
        for i in x:
            h ^= portable_hash(i)
            h *= 1000003
            h &= sys.maxsize
        h ^= len(x)
        if h == -1:
            h = -2
        return int(h)
    return hash(x)

def partitionBy(numPartitions, partitionFunc=portable_hash):
    """
    Return a copy of the RDD partitioned using the specified partitioner.

    >>> pairs = sc.parallelize([1, 2, 3, 4, 2, 4, 1]).map(lambda x: (x, x))
    >>> sets = pairs.partitionBy(2).glom().collect()
    >>> len(set(sets[0]).intersection(set(sets[1])))
    0
    """
    #if numPartitions is None:
    #    numPartitions = self._defaultReducePartitions()
    #partitioner = Partitioner(numPartitions, partitionFunc)
    #if self.partitioner == partitioner:
    #    return self

    # Transferring O(n) objects to Java is too expensive.
    # Instead, we'll form the hash buckets in Python,
    # transferring O(numPartitions) objects to Java.
    # Each object is a (splitNumber, [objects]) pair.
    # In order to avoid too huge objects, the objects are
    # grouped into chunks.
    outputSerializer = self.ctx._unbatched_serializer

    limit = (_parse_memory(self.ctx._conf.get(
        "spark.python.worker.memory", "512m")) / 2)
#======================================================================
'''


if __name__ == '__main__' :
   
   values_sep = "<>"
   keyval_sep = "::"
   b_index = "a_thaliana_t10"
   
   bowtie = prog_bowtie(values_sep, keyval_sep, b_index)
   print(bowtie.run_bowtie('CAAAGACTCATATGGACTTTGG'))

'''
http://stackoverflow.com/questions/30010939/python-subprocess-popen-error-no-such-file-or-directory
'''
