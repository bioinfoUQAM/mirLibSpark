'''
program: utils.py
author: Chao-Jung Wu
author: M.A.Remita
date: 2017-03-25
version: 1.00.01
'''

import os
import re
import subprocess
import sys

def validate_options(paramDict):
  input_type = paramDict['input_type']
  adapter = paramDict['adapter']

  if input_type == 'a' and not adapter == 'none':
    sys.stderr.write("The adapter option must be 'none' for the input_type_a.")
    sys.exit()

def makedirs_reps (reps):
  for rep in reps:
    if not os.path.exists(rep):
      os.makedirs(rep)

def find_RNAfold_path ():
  #= colosse = /software6/bioinfo/apps/mugqic_space/software/ViennaRNA/ViennaRNA-2.1.8/bin/
  proc = subprocess.Popen(['which RNAfold'], stdout=subprocess.PIPE, shell=True)
  (out, err) = proc.communicate()
  path_RNAfold = out[:-8]
  return path_RNAfold
	  
# Configure a spark context
def pyspark_configuration(appMaster, appName, masterMemory, execMemory, execCores):
  from pyspark import SparkConf, SparkContext
  myConf = SparkConf()
  #myConf.setMaster(appMaster) #= 'local[2]'
  myConf.setAppName(appName)  #= 'mirLibHadoop'
  myConf.set("spark.driver.memory", masterMemory)
  myConf.set("spark.executor.memory", execMemory) #= '4g'
  myConf.set("spark.cores.max", execCores) #= the maximum amount of CPU cores to request for the application from across the cluster (not from each machine)
  #myConf.set('spark.dynamicAllocation.enabled') # testing now --> fail on colosse
  
  #myConf.set("spark.yarn.am.memory", masterMemory) # for yarn
  #myConf.set("spark.executor.instances", execNb) # for yarn
  #myConf.set("spark.yarn.am.cores", execCores) # for yarn

  # other keys: "spark.master" = 'spark://5.6.7.8:7077'
  #             "spark.executor.memory"
  #             "spark.driver.cores"
  #             'spark.executor.cores'
  #             "spark.default.parallelism"
  return SparkContext(conf = myConf)

# Convert a file to hadoop file
# defunct
def convertTOhadoop(rfile, hdfsFile):
  print('if pre-existing in hdfs, the file would be deleted before the re-distribution of a new file with the same name.\n')
  os.system('hadoop fs -rm ' + hdfsFile) # force delete any pre-existing file in hdfs with the same name.
  os.system('hadoop fs -copyFromLocal ' + rfile + ' ' + hdfsFile)

# Convert a fasta file into a key value file
# defunct
def covert_fasta_to_KeyValue(infile, outfile):
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
      print >>fh_out, k + ':' + v
  fh_out.close()

#= Convert a seq abundance file into a key value file
# defunct
def convert_seq_freq_file_to_KeyValue(infile, outfile, v_sep):
  fh = open (infile, 'r')
  fh_out = open (outfile, 'w')

  dict_sRNA = {}
  i = 1
  
  for line in fh:
    data = line.rstrip('\n').split('\t')
    value = data[0] + v_sep + data[1]
    
    #= check if the read was treated before (redundancy)
    if data[0] not in dict_sRNA:
      dict_sRNA[data[0]] = 1
      print >>fh_out, value
      
  fh.close()
  fh_out.close()

def convert_fastq_file_to_KeyValue(infile, outfile):
  '''
  1: @SEQ_ID
  2: GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  3: +
  4: !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
  '''
  fh = open (infile, 'r')
  fh_out = open (outfile, 'w')
  i = 1
  for line in fh:
    if i == 1 or i == 3: 
      i += 1
      continue
    elif i == 2:
      seq = line.rstrip('\n')
      #print >>fh_out, seq
      i += 1
    else:
      quality = line.rstrip('\n')
      print >>fh_out, seq + '\t' + quality
      i = 1
      continue
  fh.close()
  fh_out.close()

def trim_adapter (seq, ad):
  while len(ad) > 0:
    len_ad = len(ad)
    if seq[-len_ad:] == ad:
      seq = seq[:-len_ad]
      return seq
    ad = ad[:-1]
  return seq

def getRevComp (seq):
  from string import maketrans
  
  intab = "ACGT"
  outab = "TGCA"
  trantab = maketrans(intab, outab)
  n_seq = seq.translate(trantab)
  return n_seq[::-1]

def tr_T_U (seq):
  from string import maketrans
  trantab = maketrans("T", "U")
  return seq.translate(trantab)

def tr_U_T (seq):
  from string import maketrans
  trantab = maketrans("U", "T")
  return seq.translate(trantab)

#
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

def getGenome (genome_path, file_ext):
  genome = dict()
  
  files = [each for each in os.listdir(genome_path) if each.endswith(file_ext)]
  for namefile in files :
    file = genome_path+namefile
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
    if not line.startswith("#"):
      data = line.rstrip('\r\n').split(sep)
      paramDict[data[0]] = data[1]
  
  return paramDict

def containsOnlyOneLoop (folding):
    '''
    Return True if the folding structure contains only one loop.
    A loop is definded as follows:
    When n > 0:
        loop = '(' * n1 + '.' * n2 + ')' * n3
    n1, n2, n3 don't need to be the same,
    loop1, loop2 don't need to be the same
    loop1 and loop2 are separated by a few residues composed of '.' and/or '(' and/or ')'

    @param      folding     RNAfold structure
                            ) or (  : pairing
                            .       : mismatch
    @return     True or False
    @post       folding1 = '((((((........))))))....)))'                    #True
                folding2 = '...((((((...........(((((....)))))).....))))))' #True
                folding3 = '...((.....)).....)))))..((((((...((......'      #True
                folding4 = '....((((...)))...((((......))))).....'          #False
                folding5 = '...((.....)).....))))).....((...))...'          #False
    '''
    m = re.search(r'[(]+[.]+[)]+[.()]+[(]+[.]+[)]+', folding)
    if m: return False
    return True

def writeToFile (results, outfile):
    ''' old : elem = (id, [seq, frq, nbloc, [bowtie], [pri_miRNA], [pre_miRNA]])
        new : elem = (seq, [frq, nbloc, [bowtie], [pri_miRNA], [pre_miRNA]])
    '''
    fh_out = open (outfile, 'w')    

    for elem in results :
      # ID = elem[0]#
      values = elem[1]
      miRNAseq = elem[0]
      frq = values[0]
      bowtie = values[2]
      strand = bowtie[0]#
      chromo = bowtie[1]#
      posgen = bowtie[2]#
      pre_miRNA_records = values[4]
      pre_miRNA_seq = pre_miRNA_records[0]#
      struc = pre_miRNA_records[2]#
      #mirCheck = pre_miRNA_records[3]#
      #fbstart = pre_miRNA_records[4]#
      #fbstop = pre_miRNA_records[5]#
      totalfrq = values[5]#

      #miRanda = values[6]#
      miRanda = "[TO_DO_target_genes]"
      
      #data = [miRNAseq, frq, strand, chromo, posgen, pre_miRNA_seq, struc, mirCheck, fbstart, fbstop, totalfrq, miRanda]
      data = [miRNAseq, frq, strand, chromo, posgen, pre_miRNA_seq, struc, totalfrq, miRanda]
      line = ''
      
      for d in data:
        line += str(d) + '\t'
      line = line.rstrip('\t')
      
      print >> fh_out, line
    
    fh_out.close()

def writeTimeLibToFile (timeDict, outfile, appId, paramDict):
  import datetime
  
  totalTimeSec = reduce((lambda x,y : x + y), timeDict.values())
  totalTimeHMS = str(datetime.timedelta(seconds=totalTimeSec))
  totalTimeSec = format(totalTimeSec, '.3f')
  
  fh_out = open (outfile, 'w')
  
  print >> fh_out, "# Application ID " + appId +"\n"
  print >> fh_out, "Total \t"+totalTimeHMS+ "\t" + totalTimeSec
  
  for lib in timeDict :
    timeLibSec = timeDict[lib]
    timeLibHMS = str(datetime.timedelta(seconds=timeLibSec))
    timeLibSec = format(timeLibSec, '.3f')
    
    print >> fh_out, "Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec
    print ("Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec) ##test##
  
  print >> fh_out, "\n# SPARK configuration:"
  
  for key in paramDict :
    if key.startswith("sc_"):
      print >> fh_out, "# " + key + ": " + paramDict[key]
  
  fh_out.close()

# source : https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python/23728630#23728630
def randomStrGen (n):
  import string
  import random
  
  return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(n))

def get_nonMirna_coors (infile):
  #infile = '../dbs/TAIR10_ncRNA_CDS.gff'
  idnb = 0
  d_ncRNA_CDS = {} #= {1: ['+', 'Chr5', '26939753', '26939884'], 2: ['+', 'Chr5', '26939972', '26940240'], 3: ['+', 'Chr5', '26940312', '26940396'], 4: ['+', 'Chr5', '26940532', '26940578']}
  with open (infile, 'r') as fh:
    for i in fh:
      data = i.split('\t') #= ['Chr1', 'TAIR10', 'CDS', '3760', '3913', '.', '+', '0', 'Parent=AT1G01010.1,AT1G01010.1-Protein;\n']
      chromo   = data[0]
      #genetype = data[2]   #= CDS, rRNA, snoRNA, rRNA, snoRNA, snRNA, tRNA, (gene)
      begin    = data[3]
      end      = data[4]
      strand   = data[6]
      d_ncRNA_CDS[idnb] = [chromo, strand, begin, end]
      idnb += 1
  return d_ncRNA_CDS

def get_nonMirna_list (infile, genome_path):
  # defunct
  genome = ut.getGenome (genome_path, file_ext) #= genome[chr] = sequence
  #infile = '../dbs/TAIR10_ncRNA_CDS.gff'
  l_non_miRNA = [] #= ['TGGATTTATGAAAGACGAACAACTGCGAAA']
  with open (infile, 'r') as fh:
    for i in fh:
      data = i.split('\t') #= ['Chr1', 'TAIR10', 'CDS', '3760', '3913', '.', '+', '0', 'Parent=AT1G01010.1,AT1G01010.1-Protein;\n']
      chromo   = data[0]
      if chromo == 'ChrC': chromo = 'chloroplast'
      if chromo == 'ChrM': chromo = 'mitochondria'
      begin    = int(data[3])
      end      = int(data[4])
      strand   = data[6]
      seq = genome[chromo][begin:end+1]
      if strand == '-':
        seq = ut.getRevComp (seq)
      l_non_miRNA.append(seq)
  return l_non_miRNA


def createFile_KnownNonMiRNA_from_TAIR10data (infile = 'TAIR10_GFF3_genes.gff', outfile = 'TAIR10_ncRNA_CDS.gff'):
  '''
  this function is not used in the pipeline, but users may use it to obtain their own KnonNonMiRNA from TAIR
  '''
  os.system("grep -E 'CDS|rRNA|snoRNA|snRNA|tRNA' " + infile + ' > ' + outfile)
