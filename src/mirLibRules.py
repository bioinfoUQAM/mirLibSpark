'''
program: mirLibRules.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 1.00.02

Le programme 

'''

import subprocess as sbp
import os
import utils as ut


def rearrange_rule(kv_arg, kv_sep):
  tab = kv_arg.split(kv_sep)
  return (str(tab[0]), int(tab[1]))

def rearrange_sam_rule(line):
  data = line.rstrip('\n').split('\t')
  chromo = data[2]
  posChr = int(data[3]) - 1 #= because Dr Diallo's file is not zero based
  strand = data[4]
  if strand == '1': strand = '+'
  else: strand = '-'
  seq = data[9]
  return (str(seq), [500, 1, [strand, chromo, posChr]])

def flatmap_mappings(elem) :
  '''
  ## in : ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
  ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr])
  '''
  newElems = []
  
  for mapping in elem[1][2]:
    newElem = (elem[0], [elem[1][0], elem[1][1], mapping])
    newElems.append(newElem)
  
  return newElems

class prog_dustmasker ():

  def __init__(self):
    #= The object has to be initialized in the driver program 
    #= to permit the capture of its env variables and pass them 
    #= to the subprocess in the worker nodes
    self.env = os.environ

  def dmask_filter_rule(self, elem):
    sRNAseq = str(elem[0])
    line1 = ['echo', '>seqMir\n' + sRNAseq]
    line2 = ['dustmasker']
    
    p1 = sbp.Popen(line1, stdout=sbp.PIPE, env=self.env)
    p2 = sbp.Popen(line2, stdin=p1.stdout, stdout=sbp.PIPE, env=self.env)
    p1.stdout.close()  #= Allow p1 to receive a SIGPIPE if p2 exits.
    
    output = p2.communicate()[0].rstrip('\n')
    nblines = len(output.split('\n'))
    if nblines == 1:
        return True
    return False  ##= false data will be automatically excluded in the new RDD

  def dmask_pipe_cmd(self):
    return 'dustmasker -outfmt fasta', self.env

class prog_bowtie ():

  def __init__(self, b_index):
    self.bowtie_index = b_index
    
    #= The object has to be initialized in the driver program 
    #= to permit the capture of its env variables and pass them 
    #= to the subprocess in the worker nodes
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
  
  def Bowtie_pipe_cmd (self):
    cmd = "bowtie --mm -a -v 0 --suppress 1,6,7,8 -r " + self.bowtie_index + " - "
    # cmd = "bowtie --mm -a -v 0 --suppress 1,6,7,8 " + self.bowtie_index + " - "
    
    return cmd, self.env
  
  def bowtie_rearrange_map (self, elem):
    #= u'+\tChr3\t7922370\tAAATGTAAACATCTGATCGTTTGA'
    elemTab = map(str, elem.split("\t"))
    
    #= if negative strand, get the reverse complement
    if elemTab[0] == "-" :
      elemTab[3] = ut.getRevComp(elemTab[3])
      
    return (elemTab[3], [elemTab[0], elemTab[1], int(elemTab[2])])
    
  def bowtie_freq_rearrange_rule(self, elem):
    '''
    in : ('seq',([nbLoc, [['strd','chr',posChr],..]], freq))
    out: ('seq', [freq, nbLoc, [['strd','chr',posChr],..]])
    '''
    return (elem[0], [elem[1][1], elem[1][0][0], elem[1][0][1]])

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
    
    #= all positions are zero-based
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
    
    if seq_5p == seq_3p :
      prims = [[seq_5p, pos_5p]]
    else :
      prims = [[seq_5p, pos_5p], [seq_3p, pos_3p]]
      
    return prims

  def extract_prim_rule(self, elem):
    '''
    ## in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
    ## out: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri]])
    '''
    newElems = []
    len_srna = len(elem[0])
    
    mapping = elem[1][2]
    contig = self.genome[mapping[1]]
    prims = self.extract_precursors(contig, mapping[0], mapping[2], len_srna)
      
    for prim in prims :
      newElem = (elem[0], [elem[1][0], elem[1][1], mapping, prim])
      newElems.append(newElem)
    
    return newElems

  def extract_sub_seq(self, contig, posMir, fback_start, fback_stop):
    fold_len = fback_stop - fback_start + 1
    pos = posMir - fback_start + self.pre_flank          #= 0-based
    deb = fback_start - self.pre_flank; 
    
    if deb < 0 :
      deb = 0
      pos = posMir
    
    fin = fback_stop + self.pre_flank + 1
    newSeq = contig[deb:fin]

    return [newSeq, pos]

  def extract_prem_rule(self, elem, field):
    '''
    olde : elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA, posMirPrim, Struct, mircheck, fbstart, fbstop]])
    new : elem = (seq, [frq, nbloc, [bowtie], [pri_miRNA, posMirPrim, Struct, mircheck, fbstart, fbstop]])
    '''
    priSeq = elem[1][field][0]
    posMir = int(elem[1][field][1])
    fback_start = int(elem[1][field][4])
    fback_stop = int(elem[1][field][5])
    
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
    p1.stdout.close()  #= Allow p1 to receive a SIGPIPE if p2 exits.
    
    output = p2.communicate()[0].rstrip('\n').split('\n')
    
    RNAfold_results = output[1].split(' (')
    folding = RNAfold_results[0]               #= ...(((....)))......
    # MFE = float(RNAfold_results[1][:-1])       #= minimun free energy
    
    return folding

  def RNAfold_map_rule(self, elem, field):
    '''
    old : elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA]])
    new : elem = (seq, [frq, nbloc, [bowtie], [pri_miRNA]])
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
    elem = (seq, [frq, nbloc, [bowtie], [pri_miRNA, posMirPrim, Struct]])
    '''
    len_miRNAseq = len(elem[0])
    
    pos_miRNA_start = elem[1][field][1]
    pos_miRNA_stop  = pos_miRNA_start + len_miRNAseq - 1
    folding = elem[1][field][2]
    
    mirCheck_results = self.run_mirCheck(folding, pos_miRNA_start, pos_miRNA_stop)
    
    if 'prime' in mirCheck_results[0]:
      elem[1][field].extend(mirCheck_results)
    else :
      del elem[1][field][:]
      
    return elem

class prog_dominant_profile :

  def __init__(self):
    self.env = os.environ
  
  def get_bowtie_strandchromo_dict (self, bowtie_rdd_collect):
    '''elem : (seq, [frq, nbloc, [bowties]])
    '''
    dict_bowtie_chromo_strand = {}
    
    for elem in bowtie_rdd_collect :
      bowties = elem[1][2]
      
      for bowtie in bowties :
        #= concatenate chromosome (bowtie[1]) and strand (bowtie[0])
        chromo_strand = bowtie[1] + bowtie[0]
        
        if chromo_strand not in dict_bowtie_chromo_strand.keys():
          dict_bowtie_chromo_strand[chromo_strand] = []
        
        dict_bowtie_chromo_strand[chromo_strand].append(elem)
    
    return dict_bowtie_chromo_strand
  
  def calculateTotalfrq (self, bowbloc, x, y):
    ''' old elem in bowbloc = (id, [seq, frq, nbloc, [bowties]])
        new elem in bowbloc = (seq, [frq, nbloc, [bowties]])
    '''
    totalfrq = 0
    for elem in bowbloc :
      bowties = elem[1][2]
      frq = elem[1][0]
      for bowtie in bowties :
        posgen = bowtie[2]
        if (x < posgen < y) :
          totalfrq += frq
    return totalfrq

  def profile_range (self, elem):
    ''' define x, y with pre_vld_rdd
        old : elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA], [pre_miRNA]])
        new : elem = (seq, [frq, nbloc, [bowtie], [pri_miRNA], [pre_miRNA]])
    ''' 
    posgen = elem[1][2][2]
    mirseq = elem[0]
    mirpos_on_pre = elem[1][4][1]
    preseq = elem[1][4][0]
    strand = elem[1][2][0]
    
    if strand == '+':
      x = posgen - mirpos_on_pre    #= inclusive
      y = x + len(preseq) - 1       #= inclusive
    else:
      y = posgen + len(mirseq) + mirpos_on_pre -1
      x = y-len(preseq) + 1
    return x-1, y+1                  #= exclusive  x < a < y

  def exp_profile_filter (self, elem, dict_bowtie_chromo_strand):
    ''' old : elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA], [pre_miRNA]])
        new : elem = (seq, [frq, nbloc, [bowtie], [pri_miRNA], [pre_miRNA]])
    '''
    x, y = self.profile_range (elem)
    bowtie_bloc_key = elem[1][2][1] + elem[1][2][0]  #chrom+strand
    bowbloc = dict_bowtie_chromo_strand[bowtie_bloc_key]
    totalfrq = self.calculateTotalfrq (bowbloc, x, y)
    miRNAfrq = elem[1][0]
    ratio = miRNAfrq / float(totalfrq)
    
    if ratio > 0.2 :
        return True
    return False

  def computeProfileFrq(self, elem, dict_bowtie_chromo_strand):
    x, y = self.profile_range (elem)
    bowtie_bloc_key = elem[1][2][1] + elem[1][2][0]  #=chrom+strand
    bowbloc = dict_bowtie_chromo_strand[bowtie_bloc_key]
    totalfrq = self.calculateTotalfrq (bowbloc, x, y)

    elem[1].append(totalfrq)
    return elem

class prog_miRanda ():
  def __init__ (self, Max_Score_cutoff, lower_motif_match_cutoff, upper_motif_match_cutoff, Max_Energy_cutoff, target_file, rep_tmp, miranda_exe):
    self.env = os.environ
    self.dict_seq_target = {}

    #== variables ==
    self.Max_Score_cutoff = Max_Score_cutoff
    self.lower_motif_match_cutoff = lower_motif_match_cutoff
    self.upper_motif_match_cutoff = upper_motif_match_cutoff
    self.Max_Energy_cutoff = Max_Energy_cutoff

    self.target_file = target_file
    self.rep_tmp = rep_tmp
    self.miranda_exe = miranda_exe

  def computeTargetbyMiranda (self, e):
    '''
    $miranda examples/bantam_stRNA.fasta examples/hid_UTR.fasta
    $miranda ../tmp/tmp_mirna_seq.txt ../Arabidopsis/TAIR/Genome/TAIR10_blastsets/TAIR10_cdna_20101214_updated_1cdna.fasta

    ## NOTE before disable miranda (170714): need to modify the code to use options such as -sc, -en, -go, -ge, -quiet
    '''
    if e[0] in self.dict_seq_target.keys():
      e[1].append(self.dict_seq_target[e[0]])
      return e

    tmp_file = self.rep_tmp + e[0] + 'tmp_mirna_seq_forMiranda.txt' 
    with open (tmp_file, 'w') as fh_tmp:
      print >> fh_tmp, '>x\n' + e[0]
    FNULL = open(os.devnull, 'w')
    cmd = [self.miranda_exe, tmp_file, self.target_file]

    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, stderr=FNULL, shell=False, env=self.env)
    mirandaout = sproc.communicate()[0].split('\n')
    FNULL.close()
    target_results = []
    lower_motif_match_associ, upper_motif_match_max = 0, 0

    for i in mirandaout[30:]: 
    #= because the first 30ish lines contain only program description
      if i == 'No Hits Found above Threshold': 
        ##= notTarget
        break
      if i[:2] == '>x':
        #= isTargteet_hits == ['>x', 'AT1G51370.2', '153.00', '-15.71', '2 21', '698 722', '22', '68.18%', '77.27%']
        hit_result = i.split('\t') 
        lower_motif_match_current = hit_result[7][:-1]
        upper_motif_match_current = hit_result[8][:-1]
        if float(upper_motif_match_current) > float(upper_motif_match_max):
          upper_motif_match_max = upper_motif_match_current
          lower_motif_match_associ = lower_motif_match_current
      if i[:3] == '>>x': 
        #= isTargteet == [Seq1, Seq2, Tot_Score, Tot_Energy, Max_Score, Max_Energy, Strand, Len1, Len2, Positions]
        target_result = i.split('\t') 
        Max_Score = float(target_result[4])
        Max_Energy = float(target_result[5])
        if Max_Score > self.Max_Score_cutoff and Max_Energy < self.Max_Energy_cutoff and float(upper_motif_match_max) > self.upper_motif_match_cutoff: #and float(lower_motif_match_associ) > self.lower_motif_match_cutoff:
          target_result.append(lower_motif_match_associ+'%')
          target_result.append(upper_motif_match_max+'%')
          target_results.append(target_result[1:])

    #= target_results == [[target1], [target2], ...]
    #= [['AT1G51370.2', '306.00', '-36.41', '153.00', '-20.70', '1', '23', '1118', ' 20 698', '84.21%', '89.47%']]
    #= [gene, total_score, total_energy, max_score, max_energy, strand, len_miRNA, len_gene, postions, lower_motif_match, upper_motif_match]

    self.dict_seq_target[e[0]] = target_results
    e[1].append(target_results)
    return e
  

class prog_miRdup ():
  def __init__ (self, rep_tmp, model, mirdup_jar, path_RNAfold):
    self.env = os.environ
    
    #= variable ==
    self.rep_tmp = rep_tmp
    self.model = model
    self.modelName = "".join(os.path.splitext(os.path.basename(model)))
    self.mirdup_jar = mirdup_jar
    self.path_RNAfold = path_RNAfold
    
  def run_miRdup (self, e):
    '''
    java -jar ../lib/miRdup_1.4/miRdup.jar -v ../lib/miRdup_1.4/testFiles/julie_sequencesToValidate2.txt -c ../lib/miRdup_1.4//model/Viridiplantae.model -r /usr/local/bin/
    '''

    tmp_file  = self.rep_tmp + 'mdup_'+ str(e[1]) +'.txt'
    pred_file = tmp_file+"."+self.modelName+".miRdup.txt"
    
    with open (tmp_file, 'w') as fh_tmp:
      print >> fh_tmp, 'seqx\t' + str(e[0][0]) + '\t' + str(e[0][1][4][0]) + '\t' + str(e[0][1][4][2])
    
    FNULL = open(os.devnull, 'w')

    #cmd = ['java', '-jar', '/home/cjwu/gitproject/mirLibHadoop/lib/miRdup_1.4/miRdup.jar', '-v', self.tmp_file, '-c', self.model, '-r', '/software6/bioinfo/apps/mugqic_space/software/ViennaRNA/ViennaRNA-2.1.8/bin/']
    cmd = ['java', '-jar', self.mirdup_jar, '-v', tmp_file, '-c', self.model, '-r', self.path_RNAfold]
    sproc = sbp.Popen(cmd, stdout=sbp.PIPE, stderr=FNULL, shell=False, env=self.env)
    #mirdupout = sproc.communicate()[0].split('\n')
    FNULL.close()

    mirdup_pred = ""
    mirdup_score = 0
    
    with open (pred_file, 'r') as fh_tmp :
      for line in fh_tmp :
        if line.startswith("#PR") :
          mirdup_pred = line.rstrip("\n").split("\t")[2]
        elif line.startswith("#SC") :
          mirdup_score = '%.2f' % round(float(line.rstrip("\n").split("\t")[2]), 2)

    e[0][1][4].append(mirdup_pred)
    e[0][1][4].append(mirdup_score)
    return e[0]


class prog_knownNonMiRNA ():
  def __init__ (self, non_miRNA):
    self.env = os.environ
    self.non_miRNA = non_miRNA
    
  def knFilterBySeq (self, e):
    #= defunc
    #= non_miRNA is a list of sequences
    seq = e[0]
    return not any(seq in s for s in self.non_miRNA)

  def knFilterByCoor (self, e):
    ##in : ('seq', [freq, nbLoc, ['strd','chr',posChr])
    #= non_miRNA is a list of genomic coordinations: {idnb, [strand, chromo, begin, end]}
    seq = e[0]
    bowtie = e[1][2]
    strd = bowtie[0]
    chr = bowtie[1]
    posBegin = int(bowtie[2])
    posEnd = posBegin + len(seq) - 1

    for coor in self.non_miRNA.values():
      nonmir_chromo = coor[1]
      if not chr == nonmir_chromo: continue
      nonmir_strand = coor[0]
      if not strd == nonmir_strand: continue
      nonmir_begin  = coor[2]
      if posBegin < nonmir_begin: continue
      nonmir_end    = coor[3]
      if not posEnd > nonmir_end: return False
    return True
      
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
