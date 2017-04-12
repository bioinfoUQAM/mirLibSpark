'''
program: mirLibRules.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 0.00.01

Le programme 

'''

import subprocess as sbp
import os

class filter_rules ():

  def __init__(self, vals_sep, kv_sep, limit_freq=5, limit_len=15, limit_nbLoc=5):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep
    self.limit_freq = limit_freq
    self.limit_len = limit_len
    self.limit_nbLoc = limit_nbLoc
  
  def lowfreq_filter_rule(self, kv_arg):
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    freq = value.split(self.values_sep)[1]
    if int(freq) < self.limit_freq:
        return False
    return True
  
  def shortlen_filter_rule(self, kv_arg):
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    RNAseq = value.split(self.values_sep)[0]
    if len(RNAseq) < self.limit_len:
        return False
    return True

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
    
    # The object has to be initialized in the driver program 
    # to permit the capture of its env variables and pass them 
    # to the subprocess in the worker nodes
    self.env = os.environ
    
  def dmask_filter_rule(self, kv_arg):
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    sRNAseq = value.split(self.values_sep)[0]
    line1 = ['echo', '>seq1\n' + sRNAseq]
    line2 = ['dustmasker']
    
    p1 = sbp.Popen(line1, stdout=sbp.PIPE, env=self.env)
    p2 = sbp.Popen(line2, stdin=p1.stdout, stdout=sbp.PIPE, env=self.env)
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
    
    # The object has to be initialized in the driver program 
    # to permit the capture of its env variables and pass them 
    # to the subprocess in the worker nodes
    self.env = os.environ
    
  def run_bowtie(self, seq):
    '''
    tmp_file = '/home/cloudera/workspace/miRNA_predictor/logOutput/bowtie_result_tmp.txt'
    seq = 'ATACGATCAACTAGAATGACAATT'
    #line = 'bowtie -a -v 0 --suppress 1,5,6,7,8 -c ' + self.bowtie_index + ' ' + seq + ' 1>' + tmp_file
    line = 'bowtie -a -v 0 --suppress 1,5,6,7,8 -c ' + self.bowtie_index + ' ' + seq + ' 1>/home/cloudera/workspace/miRNA_predictor/logOutput/bowtie_result_tmp.txt'
    os.system(line)
    fh = open (tmp_file, 'r')
    DATA = fh.readlines()
    fh.close()
    append_values = []
    if len(DATA) == 0:
        append_values = [['NA', 'NA', 'NA']]
    for line in DATA:
        append_value = line.rstrip('\n').split('\t')
        append_values += [append_value]
    return append_values
    '''

    #'''
    append_values = []
    FNULL = open(os.devnull, 'w')
    
    # cmd = 'bowtie --mm -a -v 0 --suppress 1,5,6,7,8 -c ' + self.bowtie_index + ' '+ seq  # shell=True
    cmd = ['bowtie', '--mm', '-a', '-v', '0', '--suppress', '1,5,6,7,8', '-c', self.bowtie_index, seq] # shell=False
    
    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, stderr=FNULL, shell=False, env=self.env)
    bsout = sproc.communicate()[0]
    bwout = bsout.decode("ascii").rstrip('\n')
    
    FNULL.close()
    
    if bwout :
      bwList = bwout.split('\n')
      for line in bwList:
        append_value = line.rstrip('\n').split('\t')
        append_values += [append_value]
      
    else :
      append_values = [['NA', 'NA', 'NA']]
    
    return append_values
    #'''
  
  def Bowtie_map_rule(self, kv_arg):
    
    keyvalue = kv_arg.split(self.keyval_sep)
    key = keyvalue[0]
    value = keyvalue[1]
    sRNAseq = value.split(self.values_sep)[0]
    
    append_value = self.run_bowtie(sRNAseq)
    
    kv2 = [key, [value] + append_value]
    return kv2



class prog_RNAfold ():
  
  def __init__(self, vals_sep, kv_sep):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep
    
    # The object has to be initialized in the driver program 
    # to permit the capture of its env variables and pass them 
    # to the subprocess in the worker nodes
    self.env = os.environ


  def run_RNAfold(self, longRNAseq):
    '''
    example line = echo GUGGAGCUCCUAUCAUUCCAAUGAAGGGUCUACCGGAAGGGUUUGUGCAGCUGCUCGUUCAUGGUUCCCACUAUCCUAUCUCCAUAGAAAACGAGGAGAGAGGCCUGUGGUUUGCAUGACCGAGGAGCCGCUUCGAUCCCUCGCUGACCGCUGUUUGGAUUGAAGGGAGCUCUGCAU | RNAfold
    task: echo and pipe the sequence to RNAfold
    this requires two subprocesses
    '''
    line1 = ['echo', longRNAseq]
    line2 = ['RNAfold']
    p1 = sbp.Popen(line1, stdout=sbp.PIPE)#, env=self.env)
    p2 = sbp.Popen(line2, stdin=p1.stdout, stdout=sbp.PIPE) #, env=self.env)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate()[0].rstrip('\n').split('\n')
    RNAfold_results = output[1].split(' (')
    folding = RNAfold_results[0] # ...(((....)))......
    MFE = float(RNAfold_results[1][:-1]) # minimun free energy
    return folding, MFE

  def RNAfold_map_rule(self, kv_arg):
    '''
    kv_arg = (id, [seq, frq, [bowtie], [pri_miRNA])
    key = kv_arg[0]
    values = kv_arg[1]
    #bowtie = values[2] #= [['+', 'ch5', '781234'], [...], ...]
    pri_miRNA = values[3] #= [[seq, pos], [seq, pos], ...] #= usually 2, or extremity 1, or special case 0
    new_pri_miRNA = [[seq, pos, folding], [seq, pos, folding], ...]
    
    pre_miRNA_example = 'GUGGAGCUCCUAUCAUUCCAAUGAAGGGUCUACCGGAAGGGUUUGUGCAGCUGCUCGUUCAUGGUUCCCACUAUCCUAUCUCCAUAGAAAACGAGGAGAGAGGCCUGUGGUUUGCAUGACCGAGGAGCCGCUUCGAUCCCUCGCUGACCGCUGUUUGGAUUGAAGGGAGCUCUGCAU'
    folding, MFE = run_RNAfold(pre_miRNA_example)
    '''
    #pri_miRNAs = kv_arg[1][3]
    for i in len(pri_miRNAs):
      seq = pri_miRNAs[i - 1][0]
      folding, MFE = run_RNAfold(seq)
      kv_arg[1][3][i - 1].append(folding)
    return kv_arg



#==================================================================
class prog_mirCheck ():
  
  def __init__(self, vals_sep, kv_sep):
    self.values_sep = vals_sep
    self.keyval_sep = kv_sep
    
    # The object has to be initialized in the driver program 
    # to permit the capture of its env variables and pass them 
    # to the subprocess in the worker nodes
    self.env = os.environ


  def run_mirCheck(self, folding, miRNA_start, miRNA_stop):
    '''
    example line = perl eval_mircheck.pl "((((((.((((((....).))))).)))))).........." 46 64 def
    '''
    cmd = ['perl', 'eval_mircheck.pl', folding, miRNA_start, miRNA_stop, 'def']
    FNULL = open(os.devnull, 'w')
    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, shell=False, stderr=FNULL)#, env=self.env)
    mirCheck_results = sproc.communicate()[0].rstrip('\n').split('\t') #= ['3prime', '1', '173']
    FNULL.close()
    return mirCheck_results

  def mirCheck_map_rule(self, kv_arg):  
    #keyvalue = kv_arg.split(self.keyval_sep)
    #key = keyvalue[0]
    #value = keyvalue[1]
    #sRNAseq = value.split(self.values_sep)[0] #= not ready

    #= artificial variables
    folding = '(((((((((((.(((.(((((...((((((...((..((((.(((.((.(((.((((.(((((....((((...(((.(((((...........))))).)))...))))....))))).)))).))).))..))).))))))..)))).))..))))).))).)))))))))))..'
    miRNA_start = '153'
    miRNA_stop = '173'
    
    append_value = mirCheck_results = self.run_mirCheck(folding, miRNA_start, miRNA_stop)
    fback = mirCheck_results[0] # True = ['3prime', '5prime']
    fback_start = mirCheck_results[1] #= '1'
    fback_stop = mirCheck_results[2] #= '173'
    print fback, fback_start, fback_stop

    #kv2 = [key, [value] + append_value] #= not ready
    #return kv2

  def mirCheck_filter_rule(self, kv_arg):
    #keyvalue = kv_arg.split(self.keyval_sep)
    #key = keyvalue[0]
    #value = keyvalue[1]
    #sRNAseq = value.split(self.values_sep)[0] #= not ready
    fback = '3prime' #= sudo affactation
    if fback == '3prime' or fback == '5prime':
        return True
    return False


#==================================================================

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
   
   #bowtie = prog_bowtie(values_sep, keyval_sep, b_index)
   #print(bowtie.run_bowtie('ATACGATCCAAGACGAGTCTCAAT'))
   mirCheck = prog_mirCheck(values_sep, keyval_sep)
   #print mirCheck.run_mirCheck('test')
   mirCheck.mirCheck_map_rule('test')
   print mirCheck.mirCheck_filter_rule('test')

'''
http://stackoverflow.com/questions/30010939/python-subprocess-popen-error-no-such-file-or-directory
'''
