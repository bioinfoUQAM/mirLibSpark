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
  return (str(tab[0]),[str(tab[1]),int(tab[2])])


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
    
    # cmd = 'bowtie --mm -a -v 0 --suppress 1,5,6,7,8 -c ' + self.bowtie_index + ' '+ seq  # shell=True

    cmd = ['bowtie', '--mm', '-a', '-v', '0', '--suppress', '1,5,6,7,8', '-c', self.bowtie_index, seq] # shell=False
    
    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, stderr=FNULL, shell=False, env=self.env)
    bsout = sproc.communicate()[0]
    bwout = bsout.rstrip('\n')
    
    FNULL.close()
    
    if bwout :
      bwList = bwout.split('\n')
      for line in bwList:
        map_value = line.rstrip('\n').split('\t')
        map_value[0] = str(map_value[0])
        map_value[1] = str(map_value[1])
        map_value[2] = int(map_value[2])
        mappings += [map_value]
      
    return mappings

  def Bowtie_map_rule(self, elem):
    sRNAseq = str(elem[1][0])
    append_value = self.run_bowtie(sRNAseq)
    elem[1].append(len(append_value))
    elem[1].append(append_value)
    return elem


class extract_precurosrs ():

  def __init__(self, genome_path, ext_left, ext_right, pre_flank):
    self.genome_path = genome_path
    self.ext_left = ext_left
    self.ext_right = ext_right
    self.pre_flank = pre_flank
    #
    self.genome = ut.getGenome(genome_path, ".fas")

  def extract_precursors (self, contig, strand, start_srna, len_srna):
    prims = []
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
    
    # seq_5p = ut.tr_T_U(seq_5p)
    # seq_3p = ut.tr_T_U(seq_3p)
    
    if seq_5p == seq_3p :
      prims = [[seq_5p, pos_5p]]
    else :
      prims = [[seq_5p, pos_5p], [seq_3p, pos_3p]]
      
    return prims

  def extract_prim_rule(self, elem):
    '''
    elem = (id, [seq, frq, nbloc, [bowties]])
    '''
    newElems = []
    primirnas = []
    len_srna = len(elem[1][0])
    
    for mapping in elem[1][3]:
      contig = self.genome[mapping[1]]
      prims = self.extract_precursors(contig, mapping[0], int(mapping[2]), len_srna)
      
      for prim in prims :
        newElem = (elem[0], [elem[1][0], elem[1][1], elem[1][2], mapping, prim])
        newElems.append(newElem)
    
    return newElems

  def extract_sub_seq(self, contig, posMir, fback_start, fback_stop):
    
    fold_len = fback_stop - fback_start + 1
    pos = posMir - fback_start + self.pre_flank          # 0-based
    deb = fback_start - self.pre_flank; 
    
    if deb < 0 :
      deb = 0
      pos = posMir
    
    fin = fback_stop + self.pre_flank + 1
    newSeq = contig[deb:fin]

    return [newSeq, pos]

  def extract_prem_rule(self, elem):
    '''
    elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA, posMirPrim, Struct, mircheck, fbstart, fbstop]])
    '''
    
    priSeq = elem[1][4][0]
    posMir = int(elem[1][4][1])
    fback_start = int(elem[1][4][4])
    fback_stop = int(elem[1][4][5])
    
    elem[1].append(self.extract_sub_seq(priSeq, posMir, fback_start, fback_stop))
    return elem


class prog_RNAfold ():

  def __init__(self):
    self.env = os.environ

  def run_RNAfold(self, seq):
    '''
    example line = echo GUGGAGCUCCUAUCAUUCC| RNAfold
    task: echo and pipe the sequence to RNAfold
    this requires two subprocesses
    '''
    line1 = ['echo', seq]
    line2 = ['RNAfold','--noPS', '--noLP']
    
    p1 = sbp.Popen(line1, stdout=sbp.PIPE, env=self.env)
    p2 = sbp.Popen(line2, stdin=p1.stdout, stdout=sbp.PIPE, env=self.env)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    
    output = p2.communicate()[0].rstrip('\n').split('\n')
    
    RNAfold_results = output[1].split(' (')
    folding = RNAfold_results[0]               # ...(((....)))......
    # MFE = float(RNAfold_results[1][:-1])       # minimun free energy
    
    return folding

  def RNAfold_map_rule(self, elem, field):
    '''
    elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA]])
    '''
    elem[1][field].append(self.run_RNAfold(elem[1][field][0]))
    return elem

class prog_mirCheck ():

  def __init__(self, param):
    self.param = param
    self.env = os.environ

  def run_mirCheck(self, folding, miRNA_start, miRNA_stop):
    '''
    example line = perl eval_mircheck.pl "((((((.((((((....).))))).)))))).........." 46 64 def
    '''
    cmd = ['perl', 'eval_mircheck.pl', folding, str(miRNA_start), str(miRNA_stop), self.param]
    FNULL = open(os.devnull, 'w')
    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, shell=False, stderr=FNULL, env=self.env)
    mirCheck_results = sproc.communicate()[0].rstrip('\n').split('\t') #= ['3prime', '1', '173']
    FNULL.close()
    return mirCheck_results

  def mirCheck_map_rule(self, elem, field):
    '''
    elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA, posMirPrim, Struct]])
    elem[1][field][0]
    '''
    len_miRNAseq = len(elem[1][0])
    
    pos_miRNA_start = int(elem[1][field][1])
    pos_miRNA_stop  = pos_miRNA_start + len_miRNAseq - 1
    folding = elem[1][field][2]
    
    mirCheck_results = self.run_mirCheck(folding, pos_miRNA_start, pos_miRNA_stop)
    
    if 'prime' in mirCheck_results[0]:
      elem[1][field].extend(mirCheck_results)
    else :
      del elem[1][field][:]
      
    # if not any(primirnas) : del primirnas[:]
    return elem

class prog_dominant_profile ():

  def __init__(self):
    self.env = os.environ

  def frq_sum_rule (self, bowelem_a, bowelem_b):
    frq_a = int(bowelem_a[1][1])
    frq_b = int(bowelem_b[1][1])
    bowelem_a[1][1] = frq_a + frq_b
    return bowelem_a #= Only the frq is useful, all other fields are meaningless.

  def filter_profile_position_rule (self, bowelem, coor):
    ''' process short RNAs with bowtie_rdd '''
    #= genome location of pre-miRNA
    #strand, chromo, x, y = '-', 'Chr3', 3366340, 3366440
    strand, chromo, x, y = coor[0],coor[1],coor[2],coor[3]
    #= return True if the genome location this short RNA is within the pre-miRNA
    nbloc = bowelem[1][2]
    for i in range(nbloc):
      if x < int(bowelem[1][3][i][2]) < y and chromo == bowelem[1][3][i][1] and strand == bowelem[1][3][i][0]:
        return True
    return False  #discard this elem

  def calculateTotalfrq (self, bowtie_rdd, coor):
    sRNAprofile = bowtie_rdd.filter(lambda bowelem: self.filter_profile_position_rule(bowelem, coor))
    totalfrq = sRNAprofile.reduce(self.frq_sum_rule)[1][1]
    return totalfrq


 
  def profile_range (elem):
    ''' define x, y with pre_vld_rdd'''
    posgen = elem[1][3][2] 		# already int
    mirseq = elem[1][0]
    mirpos_on_pre = elem[1][5][1] 	# already int
    preseq = elem[1][5][0]
    strand = elem[1][3][0]
    if strand == '+':
      x = posgen - mirpos_on_pre	# inclusive
      y = x + len(preseq) - 1		# inclusive
    else:
      y = posgen + len(mirseq) + mirpos_on_pre -1
      x = y-len(preseq) + 1
    return x-1, y+1	# exclusive  x < a < y


  def functionX (self, elem, bowtie_rdd):
    x, y = profile_range (elem)
    totalfrq = calculateTotalfrq (bowtie_rdd, coor)
    miRNAfrq = elem[1][1]
    ratio = miRNAfrq / totalfrq
    if ratio == 1:
        return True
    
    return True # if miRNAfrq_percentage > 0.2

    





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
