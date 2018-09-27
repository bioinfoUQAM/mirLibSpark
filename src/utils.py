'''
program: utils.py
author: Chao-Jung Wu
author: M.A.Remita
date: 2017-03-25
update: 2018-09-18
version: 1.00.01
'''
from __future__ import print_function
import os
import re
import subprocess
import sys


def validate_options(paramDict):
  '''
  wip: examine the combination of the parameters to see if it is compatible logically
  '''
  input_type = paramDict['input_type']
  adapter = paramDict['adapter']

  if input_type == 'a' and not adapter == 'none':
    sys.stderr.write("The adapter option must be 'none' for the input_type_a.")
    sys.exit()

def transpose_txt(infile, outfile):
    with open(infile, 'r') as f:
        lis = [x.rstrip('\n').split('\t') for x in f]
    fho = open (outfile, 'w')
    for x in zip(*lis):
        for y in x:
            print(y+'\t', end='', file=fho)
            #print >>fho, y+'\t', # the comma in the final signals not to print a new line in python2.x
        print('', file=fho)
        #print >>fho, ''
		
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
  myConf.setMaster(appMaster) #= 'local[2] or local[*]'
  myConf.setAppName(appName)  #= 'mirLibHadoop'
  myConf.set("spark.driver.memory", masterMemory)
  myConf.set("spark.executor.memory", execMemory) #= '4g'
  myConf.set("spark.cores.max", execCores) 
  
  # other keys: "spark.master" = 'spark://5.6.7.8:7077'
  #             "spark.driver.cores"
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
      #print >>fh_out, k + ':' + v
      print (k + ':' + v, file=fh_out)
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
      #print >>fh_out, value
      print (value, file=fh_out)
      
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
      #print >>fh_out, seq + '\t' + quality
      print (seq + '\t' + quality, file=fh_out)
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
  intab = "ACGT"
  outab = "TGCA"
  #= v1 ================================
  from string import maketrans
  trantab = maketrans(intab, outab)
  #= v2 ================================
  #trantab = str.maketrans(intab, outab)
  #=====================================
  n_seq = seq.translate(trantab)
  return n_seq[::-1]
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

def getGenome_ (genome_path, file_ext):
  ''' defunct '''
  genome = dict()
  
  files = [each for each in os.listdir(genome_path) if each.endswith(file_ext)]
  for namefile in files :
    file = genome_path+namefile
    chr = getChromosomeName(file)
    sequence = getFastaSeq(file)
    genome[chr] = sequence
    
  return genome

def getGenome (genome_path, file_ext, chromosomeName):
  genome = dict()
  
  if chromosomeName == 'All':
    files = [each for each in os.listdir(genome_path) if each.endswith(file_ext)]
  else:
    files = [ chromosomeName + file_ext]

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
    if line.startswith('message'): 
      msg = line.rstrip('\r\n')[8:].split('\\n')
      for i in msg: print(i)
    elif not line.startswith("#"):
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
    ## result in: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore'], totalfrq]) 
    fh_out = open (outfile, 'w')

    ## result out:
    line = '\t'.join('miRNAseq, frq, nbLoc, strand, chromo, posChr, mkPred, mkStart, mkStop, preSeq, posMirPre, newfbstart, newfbstop, preFold, mpPred, mpScore, totalFrq'.split(', ')) 
    #print >> fh_out, line
    print(line, file=fh_out)

    seen = []
    for elem in results :
      miRNAseq = elem[0]
      frq = elem[1][0]
      #= bowtie_result
      nbLoc = elem[1][1]
      strand = elem[1][2][0]
      chromo = elem[1][2][1]
      posChr = elem[1][2][2]
      #= pri-miRNA
      priSeq = elem[1][3][0] 
      posMirPri = elem[1][3][1] 
      priFold = elem[1][3][2] 
      #= mircheck_result
      mkPred = elem[1][3][3]
      mkStart = elem[1][3][4]
      mkStop = elem[1][3][5]  
      #= mirdup_result
      preSeq = elem[1][4][0]
      posMirPre = elem[1][4][1]
      newfbstart = int(posMirPre) + int(posMirPri) - int(mkStart) 
      newfbstop  = int(posMirPre) + int(mkStop) - int(posMirPri) 
      preFold = elem[1][4][2]
      mpPred = elem[1][4][3]
      mpScore = elem[1][4][4]
      #= premirna_range_total_small_rna_freq
      totalFrq =  elem[1][5]
      
      seeing = ''.join([str(x) for x in [miRNAseq, strand, chromo, posChr, mkPred, len(preSeq), posMirPre, newfbstart, newfbstop, mpPred, mpScore]])
      if seeing in seen: continue
      else: seen.append(seeing)
      #data = [miRNAseq, frq, nbLoc, strand, chromo, posChr, mkPred, preSeq, posMirPre, preFold, mpPred, mpScore, totalFrq]
      data = [miRNAseq, frq, nbLoc, strand, chromo, posChr, mkPred, mkStart, mkStop, preSeq, posMirPre, newfbstart, newfbstop, preFold, mpPred, mpScore, totalFrq] ##update

      line = '\t'.join([str(d) for d in data])
      #print >> fh_out, line
      print(line, file=fh_out)
    fh_out.close()


def writeTimeLibToFile (timeDict, outfile, appId, paramDict):
  import datetime
  
  totalTimeSec = reduce((lambda x,y : x + y), timeDict.values())
  totalTimeHMS = str(datetime.timedelta(seconds=totalTimeSec))
  totalTimeSec = format(totalTimeSec, '.3f')
  
  fh_out = open (outfile, 'w')
  
  #print >> fh_out, datetime.datetime.now()
  print(datetime.datetime.now(), file=fh_out)

  #print >> fh_out, "# Application ID " + appId +"\n"
  print("# Application ID " + appId +"\n", file=fh_out)

  #print >> fh_out, "Total \t"+totalTimeHMS+ "\t" + totalTimeSec
  print("Total \t"+totalTimeHMS+ "\t" + totalTimeSec, file=fh_out)


  for lib in timeDict :
    timeLibSec = timeDict[lib]
    timeLibHMS = str(datetime.timedelta(seconds=timeLibSec))
    timeLibSec = format(timeLibSec, '.3f')
    
    #print >> fh_out, "Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec
    print("Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec, file=fh_out)
    print ("Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec) #= stdout
  
  #print >> fh_out, "\n# SPARK configuration:"
  print("\n# SPARK configuration:", file=fh_out)

  for key in paramDict :
    if key.startswith("sc_"):
      #print >> fh_out, "# " + key + ": " + paramDict[key]
      print("# " + key + ": " + paramDict[key], file=fh_out)


  #print >> fh_out, "\n# MirLibSpark configuration:"
  print("\n# MirLibSpark configuration:", file=fh_out)

  for key in sorted(paramDict.keys()) :
    if not key.startswith("sc_"):
      #print >> fh_out, "# " + key + ": " + paramDict[key]
      print("# " + key + ": " + paramDict[key], file=fh_out)


  
  fh_out.close()


###########################
## WORK IN PROGRESS
###########################
def writeSummaryExpressionToFile (infiles, rep_output, appId):
  '''
  ## in : [miRNAseq, frq, nbLoc, strand, chromo, posChr, mkPred, mkStart, mkStop, preSeq, posMirPre, newfbstart, newfbstop, preFold, mpPred, mpScore, totalFrq]
  '''
  keyword = '_miRNAprediction_'
  #= .trs is a temporary extension, such file will be transposed at the end of this function
  outfile = rep_output + appId + '_summaryFreq.trs' 
  outfile2 = rep_output + appId + '_summaryBinary.trs'
  #outfile3 = rep_output + appId + '_summaryGenoLoci.txt' #= outfile3 is not in a good format yet

  fh_out = open (outfile, 'w')
  fh_out2 = open (outfile2, 'w')
  #fh_out3 = open (outfile3, 'w')
  
  #line_seen = []
  #for f in sorted(infiles):
  #  with open (rep_output + f, 'r') as fh:
  #    for line in fh:
  #      line=line.rstrip('\n')
  #      if line not in line_seen:
  #        line_seen.append(line)
  #      print >> fh_out3, line
  #fh_out3.close()

  master_predicted_distinctMiRNAs = []
  tmp_master_distinctPrecursor_infos = {}
  for f in sorted(infiles):
    with open (rep_output + f, 'r') as fh:
      ###======================
      DATA = list(set(fh.readlines()[1:]))
      ###======================
      for line in DATA:
        data = line.rstrip('\n').split('\t')
        ########################################## 
        miRNAseq = data[0]
        frq = data[1]
        nbLoc = data[2]
        strand = data[3]
        chromo = data[4]
        posChr = data[5]
        mkPred = data[6]
        mkStart = data[7]
        mkStop = data[8]  
        preSeq = data[9]
        posMirPre = data[10]
        newfbstart = data[11]
        newfbstop  = data[12]
        preFold = data[13]
        mpPred = data[14]
        mpScore = data[15]
        totalFrq =  data[16]
        #################################### 
        if miRNAseq not in master_predicted_distinctMiRNAs:
          master_predicted_distinctMiRNAs.append(miRNAseq)
        #################################### 
        key = miRNAseq + ':'+ newfbstart + newfbstop
        if key not in tmp_master_distinctPrecursor_infos.keys():
          infos = [miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore]
          tmp_master_distinctPrecursor_infos[key] = infos
        #################################### 

  dictLibSeqFreq = {}
  for f in sorted(infiles):
    libname = f.split(keyword)[1][:-4]
    dictLibSeqFreq[libname] = []
    tmpDict = {}
    with open (rep_output + f, 'r') as fh:
      fh.readline()
      for line in fh:
        data = line.rstrip('\n').split('\t')
        miRNAseq = data[0]
        freq = int(data[1])
        tmpDict[miRNAseq] = freq
    for e in master_predicted_distinctMiRNAs:
      if e in tmpDict.keys(): dictLibSeqFreq[libname].append(tmpDict[e])
      else: dictLibSeqFreq[libname].append(0)

  seqListLine = 'miRNA\t' + '\t'.join(master_predicted_distinctMiRNAs)

  #print >> fh_out, seqListLine
  print(seqListLine, file=fh_out)

  #print >> fh_out2, seqListLine
  print(seqListLine, file=fh_out2)

  for k in sorted(dictLibSeqFreq.keys()):
    v = dictLibSeqFreq[k]
    line = k + '\t'
    for i in v: line += str(i) + '\t'
    line = line.rstrip('\t')
    #print >> fh_out, line
    print(line, file=fh_out)


  for k in sorted(dictLibSeqFreq.keys()):
    v = dictLibSeqFreq[k]
    line = k + '\t'
    for i in v:
      if i > 0: i = 1
      line += str(i) + '\t'
    line = line.rstrip('\t')
    #print >> fh_out2, line
    print(line, file=fh_out2)


  master_distinctPrecursor_infos = []
  for k in sorted(tmp_master_distinctPrecursor_infos.keys()):
    v = tmp_master_distinctPrecursor_infos[k]
    master_distinctPrecursor_infos.append(v)

  fh_out.close()
  fh_out2.close()

  import os
  infile = outfile
  outfile = infile[:-4] + '.txt'
  transpose_txt(infile, outfile)
  os.remove(infile) 

  infile = outfile2
  outfile = infile[:-4] + '.txt'
  transpose_txt(infile, outfile)
  os.remove(infile) 

  #return master_predicted_distinctMiRNAs #sorted(master_distinctMiRNAs_infos)
  return sorted(master_distinctPrecursor_infos)


def writeTargetsToFile (mirna_and_targets, rep_output, appId):
  '''
  targets = [miRNAseq, [targets], mirnazipindex]
  '''
  outfile = rep_output + appId + '_mirna_and_targets.txt'
  fh_out = open (outfile, 'w')
  outfile2 = rep_output + appId + '_mirna_and_topscoredTargets.txt'
  fh_out2 = open (outfile2, 'w')
  topKscored = 5

  for i in mirna_and_targets:
    print >> fh_out, i[0], i[1], str(i[2]).zfill(4)

    mirnaseq = i[0]
    mirna_uindex = i[2]
    targetcollect = []
    seen = []
    score = 1000000000000
    count = 0
    for t in i[1]:
      target = t[0].split('.')[0]
      score_cur = t[1]
      if score_cur < score:
        count += 1
        score = score_cur
      if count < (topKscored +1) and target not in seen:
        targetcollect.append( target + ' ('+ str(score_cur) +')' )
        seen.append(target)
    #data = [mirnaseq, targetcollect[0], ','.join(targetcollect)]
    data = [mirnaseq, ','.join(targetcollect)]
    line = '\t'.join(data)
    #print >> fh_out2, line
    print(line, file=fh_out2)


  fh_out.close()
  fh_out2.close()

def __charge_gene_vs_pathway_file (infile):
  d_gene_pathway = {} # {'gene': 'p1,p2,p3'}
  with open (infile) as fh: DATA = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  for i in DATA:
    gene = i[0]
    p = i[1]
    if gene not in d_gene_pathway.keys(): d_gene_pathway[gene] = p
    else: d_gene_pathway[gene] += ',' + p
  return d_gene_pathway

def annotate_target_genes_with_KEGGpathway (gene_vs_pathway_file, rep_output, appId):
  d_gene_pathway = __charge_gene_vs_pathway_file (gene_vs_pathway_file)
  infile = rep_output + appId + '_mirna_and_topscoredTargets.txt'
  outfile = rep_output + appId + '_mirna_and_topscoredTargetsKEGGpathway.txt'
  fh_out = open (outfile, 'w')
  newDATA = []
  with open (infile) as fh: DATA = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  for i in DATA:
    mirna = i[0]
    data = [mirna]
    for method in i[1:]: #=method = top1, top5scored
      genes = method.split(',')
      geneAnnotation = []
      for g in genes:
        g = g.split(' (')[0].split('.')[0]
        if g in d_gene_pathway.keys(): 
          p = d_gene_pathway[g]
          geneAnnotation.append(p)
      data.append( ','.join(geneAnnotation) )
    newDATA.append(data)
  for i in newDATA:
    line = '\t'.join(i)
    #print >> fh_out, line
    print(line, file=fh_out)

  fh_out.close()





#############################
#############################
def write_index (data, rep_output, appId):
  #=serial, miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore
  outfile = rep_output + appId + '_precursorindex.txt'
  fh_out = open (outfile, 'w')
  for i in data:
    i[0] = str(i[0]).zfill(4)
    line = '\t'.join( [str(x) for x in i] )
    #print >> fh_out, line
    print(line, file=fh_out)

  fh_out.close()
  write_html (data, rep_output, appId)


def write_html (DATA, rep_output, appId):
  #=serial, miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore
  infile = rep_output + appId + '_precursorindex.txt'
  outfile = rep_output + appId +'_precursorindex.html'
  fh_out=open(outfile,'w')
  '''
  l='<html>\n<head>\n<style>';print >> fh_out, l
  l='table{\nfont-family:arial,sans-serif;\nborder-collapse:collapse;\nwidth:90%;\nmargin:auto;\n}';print >> fh_out, l
  l='td,th{\nborder:1pxsolid#dddddd;\ntext-align:left;\npadding:8px;\nmax-width:200px;\nword-break:break-all;\n}';print >> fh_out, l
  l='tr:nth-child(even){background-color:#dddddd;}';print >> fh_out, l
  l='img.heightSet{max-width:100px;max-height:200px;}';print >> fh_out, l
  l='img:hover{box-shadow:002px1pxrgba(0,140,186,0.5);}';print >> fh_out, l
  l='</style></head><body>';print >> fh_out, l
  l='<div GenomicPre><h2>miRNAs and their genomic precursors</h2><table>';print >> fh_out, l
  l='  <tr>';print >> fh_out, l
  l='    <th>Serial</th>';print >> fh_out, l
  #l='    <th>NewID</th>';print >> fh_out, l
  l='    <th>image</th>';print >> fh_out, l
  l='    <th>miRNA.pre-miRNA.Structure</th>';print >> fh_out, l
  l='    <th>Coordination</th>';print >> fh_out, l
  l='  </tr>';print >> fh_out, l
  '''

  l='<html>\n<head>\n<style>';print(l, file=fh_out)
  l='table{\nfont-family:arial,sans-serif;\nborder-collapse:collapse;\nwidth:90%;\nmargin:auto;\n}';print(l, file=fh_out)
  l='td,th{\nborder:1pxsolid#dddddd;\ntext-align:left;\npadding:8px;\nmax-width:200px;\nword-break:break-all;\n}';print(l, file=fh_out)
  l='tr:nth-child(even){background-color:#dddddd;}';print(l, file=fh_out)
  l='img.heightSet{max-width:100px;max-height:200px;}';print(l, file=fh_out)
  l='img:hover{box-shadow:002px1pxrgba(0,140,186,0.5);}';print(l, file=fh_out)
  l='</style></head><body>';print(l, file=fh_out)
  l='<div GenomicPre><h2>miRNAs and their genomic precursors</h2><table>';print(l, file=fh_out)
  l='  <tr>';print(l, file=fh_out)
  l='    <th>Serial</th>';print(l, file=fh_out)
  #l='    <th>NewID</th>';print(l, file=fh_out)
  l='    <th>miRNA.pre-miRNA.Structure</th>';print(l, file=fh_out)
  l='    <th>Coordination</th>';print(l, file=fh_out)
  l='  </tr>';print(l, file=fh_out)


  ## loop start
  for i in DATA:
    i = [str(x) for x in i]
    serial = i[0]
    newid = 'na'	
    mirna = i[1]	
    chromo = i[3]	
    poschromo = i[4]	
    preseq = i[5]	
    #pospre = i[6]	
    structure = i[7]	
    #mircheck = i[8:11]	
    strand = i[2]
 
    path = rep_output + appId + '_' + serial.zfill(4) + '_' + chromo + '_' + poschromo + '.jpg'
    '''
    l='  <tr>';print >> fh_out, l
    l="    <td rowspan=3 style='width: 120px;'><strong>"+ serial + "</strong></td>";print >> fh_out, l
    #l="    <td rowspan=3 style='width: 120px;'>"+newid+"</td>";print >> fh_out, l
    l="    <td rowspan=3 > <a href="+ path + " target='_blank'><img class='heightSet' src="+ path + " alt='Struture'></a></td>";print >> fh_out, l
    l="    <td style='width: 400px;'>"+ mirna + "</td>";print >> fh_out, l
    l="    <td>"+ chromo + ":"+ poschromo + " ["+ strand + "] </td>";print >> fh_out, l
    l="  </tr>";print >> fh_out, l
    l="  <tr>";print >> fh_out, l
    l="    <td colspan=2 style='font-family:monospace'>"+ preseq + "</td>";print >> fh_out, l
    l="  </tr>";print >> fh_out, l
    l="  <tr>";print >> fh_out, l
    l="    <td colspan=2 style='font-family:monospace'>"+ structure + "</td>";print >> fh_out, l
    l="  </tr>";print >> fh_out, l
  ## loop end
  
  l="</table></div>";print >> fh_out, l
  l="</body></html>";print >> fh_out, l
  '''


    l='  <tr>';print(l, file=fh_out)
    l="    <td rowspan=3 style='width: 120px;'><strong>"+ serial + "</strong></td>";print(l, file=fh_out)
    #l="    <td rowspan=3 style='width: 120px;'>"+newid+"</td>";print(l, file=fh_out)
    l="    <td rowspan=3 > <a href="+ path + " target='_blank'><img class='heightSet' src="+ path + " alt='Struture'></a></td>";print(l, file=fh_out)
    l="    <td style='width: 400px;'>"+ mirna + "</td>";print(l, file=fh_out)
    l="    <td>"+ chromo + ":"+ poschromo + " ["+ strand + "] </td>";print(l, file=fh_out)
    l="  </tr>";print(l, file=fh_out)
    l="  <tr>";print(l, file=fh_out)
    l="    <td colspan=2 style='font-family:monospace'>"+ preseq + "</td>";print(l, file=fh_out)
    l="  </tr>";print(l, file=fh_out)
    l="  <tr>";print(l, file=fh_out)
    l="    <td colspan=2 style='font-family:monospace'>"+ structure + "</td>";print(l, file=fh_out)
    l="  </tr>";print(l, file=fh_out)
  ## loop end
  
  l="</table></div>";print(l, file=fh_out)
  l="</body></html>";print(l, file=fh_out)


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
  genome = ut.getGenome (genome_path, file_ext, 'All') #= genome[chr] = sequence
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


