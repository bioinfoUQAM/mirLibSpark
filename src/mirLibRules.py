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
import utils as ut


def rearrange_rule(kv_arg, kv_sep):
  tab = kv_arg.split(kv_sep)
  return (tab[0],[tab[1],tab[2]])

class prog_dustmasker ():
  
  def __init__(self):
    # The object has to be initialized in the driver program 
    # to permit the capture of its env variables and pass them 
    # to the subprocess in the worker nodes
    self.env = os.environ
    
  def dmask_filter_rule(self, elem):
    sRNAseq = str(elem[1][0])
    line1 = ['echo', '>seqMir\n' + sRNAseq]
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
  
  def __init__(self, b_index):
    self.bowtie_index = b_index
    
    # The object has to be initialized in the driver program 
    # to permit the capture of its env variables and pass them 
    # to the subprocess in the worker nodes
    self.env = os.environ
    
  def run_bowtie(self, seq):
    mappings = []
    FNULL = open(os.devnull, 'w')
    
    # self.cmd = 'bowtie --mm -a -v 0 --suppress 1,5,6,7,8 -c ' + self.bowtie_index + ' '+ seq  # shell=True

    cmd = ['bowtie', '--mm', '-a', '-v', '0', '--suppress', '1,5,6,7,8', '-c', self.bowtie_index, seq] # shell=False
    
    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, stderr=FNULL, shell=False, env=self.env)
    bsout = sproc.communicate()[0]
    bwout = bsout.decode("ascii").rstrip('\n')
    
    FNULL.close()
    
    if bwout :
      bwList = bwout.split('\n')
      for line in bwList:
        map_value = line.rstrip('\n').split('\t')
        mappings += [map_value]
      
    return mappings

  
  def Bowtie_map_rule(self, elem):
    sRNAseq = str(elem[1][0])
    append_value = self.run_bowtie(sRNAseq)
    elem[1].append(append_value)
    return elem


#== sudo code ==
class sudo ():
  def __init__(self):
    self.env = os.environ

  def sudo_long (self, kv_arg):
    '''
    kv_arg = (id, [seq, frq, [bowtie], [pri_miRNA])
    kv_arg_after = (id, [seq, frq, [bowtie], [pri_miRNA])
    '''
    pre_miRNA_example = 'GUGGAGCUCCUAUCAUUCCAAUGAAGGGUCUACCGGAAGGGUUUGUGCAGCUGCUCGUUCAUGGUUCCCACUAUCCUAUCUCCAUAGAAAACGAGGAGAGAGGCCUGUGGUUUGCAUGACCGAGGAGCCGCUUCGAUCCCUCGCUGACCGCUGUUUGGAUUGAAGGGAGCUCUGCAU'
    #bowtie = kv_arg[1][2]
    #pri_miRNA = kv_arg[1][3]
    kv_arg[1].append([])
    nbLoc = len(kv_arg[1][2])
    for i in range(nbLoc):
      kv_arg[1][3][i].append( [pre_miRNA_example, '153'] )
    return kv_arg  
#===============

class prog_RNAfold ():
  
  def __init__(self):
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
    folding = RNAfold_results[0]               # ...(((....)))......
    MFE = float(RNAfold_results[1][:-1])       # minimun free energy
    return folding, MFE

  def RNAfold_map_rule(self, kv_arg):
    '''
    kv_arg = (id, [seq, frq, [bowtie], [pri_miRNA])
    pre_miRNA_example = 'GUGGAGCUCCUAUCAUUCCAAUGAAGGGUCUACCGGAAGGGUUUGUGCAGCUGCUCGUUCAUGGUUCCCACUAUCCUAUCUCCAUAGAAAACGAGGAGAGAGGCCUGUGGUUUGCAUGACCGAGGAGCCGCUUCGAUCCCUCGCUGACCGCUGUUUGGAUUGAAGGGAGCUCUGCAU'
    folding, MFE = run_RNAfold(pre_miRNA_example)
    '''
    nbLoc = len(kv_arg[1][2])
    pri_miRNAs = kv_arg[1][3]
    for i in range(nbLoc):
      nb_pri_miRNA_OF_thisLoc = len(kv_arg[1][3][i])
      for j in range(nb_pri_miRNA_OF_thisLoc):
        pri_seq = kv_arg[1][3][i][0]
        folding, MFE = self.run_RNAfold(pri_seq)
        kv_arg[1][3][i].append(folding)
    return kv_arg
#==================================================================
class prog_mirCheck ():
  
  def __init__(self):    
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
    '''
    kv_arg = (id, [seq, frq, [bowtie], [pri_miRNA])
    key = kv_arg[0]
    values = kv_arg[1]
    #bowtie = values[2] #= [['+', 'ch5', '781234'], [...], ...]
    pri_miRNA = values[3] #= [[seq, pos], [seq, pos], ...] #= usually 2, or extremity 1, or special case 0
    new_pri_miRNA = [[seq, pos, folding], [seq, pos, folding], ...]
    '''
    kv_arg[1].append([])    #= pre_miRNAs = kv_arg[1][4]
    len_miRNAseq = len(kv_arg[1][0])
    nbLoc = len(kv_arg[1][2])
    pri_miRNAs = kv_arg[1][3]
    for i in range(nbLoc):
      nb_pri_miRNA_OF_thisLoc = len(kv_arg[1][3][i]) #########################
      for j in range(nb_pri_miRNA_OF_thisLoc):
        pri_seq = kv_arg[1][3][i][0]
        pos_miRNA_start = kv_arg[1][3][i][1]
        folding = kv_arg[1][3][i][2]
        pos_miRNA_stop = str(  int(pos_miRNA_start) + len_miRNAseq - 1  )

        mirCheck_results = self.run_mirCheck(folding, pos_miRNA_start, pos_miRNA_stop)
        #fback = mirCheck_results[0] # True = ['3prime', '5prime']
        #fback_start = mirCheck_results[1] #= '1'
        #fback_stop = mirCheck_results[2] #= '173'
        #print fback, fback_start, fback_stop
        ###kv_arg[1][4][i].append(mirCheck_results)
        #kv_arg[1][4][i] += mirCheck_results
        
    return kv_arg


    #= artificial variables
    folding = '(((((((((((.(((.(((((...((((((...((..((((.(((.((.(((.((((.(((((....((((...(((.(((((...........))))).)))...))))....))))).)))).))).))..))).))))))..)))).))..))))).))).)))))))))))..'
    miRNA_start = '153'
    miRNA_stop = '173'
    
    append_value = mirCheck_results = self.run_mirCheck(folding, miRNA_start, miRNA_stop)
    fback = mirCheck_results[0] # True = ['3prime', '5prime']
    fback_start = mirCheck_results[1] #= '1'
    fback_stop = mirCheck_results[2] #= '173'
    print fback, fback_start, fback_stop


  def mirCheck_filter_rule(self, kv_arg):
    #keyvalue = kv_arg.split(self.keyval_sep)
    #key = keyvalue[0]
    #value = keyvalue[1]
    #sRNAseq = value.split(self.values_sep)[0] #= not ready
    #fback = '3prime' #= sudo affactation
    #if fback == '3prime' or fback == '5prime':
    #    return True
    #return False
    pre_miRNAs = kv_arg[1][4] #= [['3prime', '1', '173'], [...], ...]
    nbLoc = len(kv_arg[1][2])
    for i in range(nbLoc):
      nb_pre_miRNA_OF_thisLoc = len(kv_arg[1][4][i])
      for j in range(nb_pre_miRNA_OF_thisLoc):
        fback = kv_arg[1][4][i][j][0]
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

class extract_precurosrs ():
  
  def __init__(self, genome_path, ext_left, ext_right):
    self.genome_path = genome_path
    self.ext_left = ext_left
    self.ext_right = ext_right
    #
    self.genome = 0 ###
    
  #def getGenome ():
    #return 0
    
  def extract_precursor (self, contig, strand, start_srna, len_srna):
    ext_left = self.ext_left
    ext_right = self.ext_right
    
    # all positions are zero-based
    s1 = start_srna - ext_left
    pos_5p = ext_left
    if s1 < 0:
      s1 = 0
      pos_5p = start_srna
      ext_left = start_srna
    
    t_len1 = len_srna + ext_left + ext_right
    
    s2 = start_srna - ext_right
    pos_3p = ext_right
    if s2 < 0:
      s2 = 0
      pos_3p = start_srna
      ext_right = start_srna
      
    t_len2 = len_srna + ext_left + ext_right
    
    seq_5p = contig[s1:(s1+t_len1)]
    seq_3p = contig[s2:(s2+t_len2)]
    
    if strand == "-":
      seq_5p = ut.getRevComp(seq_5p);
      pos_5p = len(seq_5p) - pos_5p - len_srna
      
      seq_3p = ut.getRevComp(seq_3p);
      pos_3p = len(seq_3p) - pos_3p - len_srna
    
    seq_5p = ut.tr_T_U(seq_5p)
    seq_3p = ut.tr_T_U(seq_3p)
    
    return [[seq_5p, pos_5p], [seq_3p, pos_3p]]
    
  def extract_prec_rule(self, elem):
    # (u'000000001', [u'GAGGTTGGACAAGGCTTT', u'549', [[u'+', u'Chr2', u'1647574']]])
    primirnas = []
    len_srna = len(elem[1][0])
    
    for mapping in elem[1][2]:
      contig = self.genome[mapping[1]]
      primirnas.append(extract_precursor(contig, mapping[0], mapping[2], len_srna))
    
    elem[1].append(primirnas)
    return elem
  
  
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
