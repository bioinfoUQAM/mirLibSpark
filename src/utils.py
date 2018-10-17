'''
program: utils.py
author: Chao-Jung Wu
author: M.A.Remita
date: 2017-03-25
update: 2018-09-18
version: 1.00.01

support python3 print syntax
'''
from __future__ import print_function
import os
import re
import subprocess
import sys

import os.path
from os import listdir



def validate_options(paramDict):
  '''
  wip: examine the combination of the parameters to see if it is compatible logically
  '''
  input_type = paramDict['input_type']
  adapter = paramDict['adapter']
  rep_input = paramDict['input_path']
  diffguide_file = paramDict['diffguide_file']
  perform_differnatial_analysis = paramDict['perform_differnatial_analysis']
  perform_KEGGpathways_enrichment_analysis= paramDict['perform_KEGGpathways_enrichment_analysis']



  if input_type == 'a' and not adapter == 'none':
    sys.stderr.write("The adapter option must be 'none' for the input_type_a.\nExit the program.")
    sys.exit()

  if perform_KEGGpathways_enrichment_analysis == 'yes' and perform_differnatial_analysis == 'no':
    sys.stderr.write("KEGG pathway enrichment analysis must be done after differential expression analysis.\nExit the program.")
    sys.exit()



  #= verify if input folder contain all files requisted by diffguide file
  if perform_differnatial_analysis == 'yes':
    infiles = [f.split('.')[0] for f in listdir(rep_input) if os.path.isfile(os.path.join(rep_input, f))]
    #testInfiles = [f.split('.')[0] for f in infiles]
    diffguide, neededInfiles = __read_diffguide(diffguide_file)
    for infile in neededInfiles:
      if infile not in infiles: 
        sys.stderr.write('One or more input files requested by diffguide_file are not provided in the input folder.\nExit the program.')
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
  myConf.setAppName(appName)  #= 'mirLibSpark'
  myConf.set("spark.driver.memory", masterMemory)
  myConf.set("spark.executor.memory", execMemory) 
  myConf.set("spark.cores.max", execCores) 
  
  # other keys: "spark.master" = 'spark://5.6.7.8:7077'
  #             "spark.driver.cores"
  #             "spark.default.parallelism"
  return SparkContext(conf = myConf)

def convertTOhadoop(rfile, hdfsFile):
  '''
  Convert a file to hadoop file
  defunct
  '''
  print('if pre-existing in hdfs, the file would be deleted before the re-distribution of a new file with the same name.\n')
  os.system('hadoop fs -rm ' + hdfsFile) # force delete any pre-existing file in hdfs with the same name.
  os.system('hadoop fs -copyFromLocal ' + rfile + ' ' + hdfsFile)

def covert_fasta_to_KeyValue(infile, outfile):
  '''
  Convert a fasta file into a key value file
  defunct
  '''
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
      print (k + ':' + v, file=fh_out)
  fh_out.close()

def convert_seq_freq_file_to_KeyValue(infile, outfile, v_sep):
  '''
  Convert a seq abundance file into a key value file
  defunct
  '''
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
      i += 1
    else:
      quality = line.rstrip('\n')
      print (seq + '\t' + quality, file=fh_out)
      i = 1
      continue
  fh.close()
  fh_out.close()

def find_str(s, char):
    ''' zero based indexing '''
    index = 0
    if char in s:
        c = char[0]
        for ch in s:
            if ch == c:
                if s[index:index+len(char)] == char:
                    return index
            index += 1
    return -1000

'''
#ORIGINAL FUNCTION
'''
def trim_adapter__ (seq, ad):
  while len(ad) > 0:
    len_ad = len(ad)
    if seq[-len_ad:] == ad:
      seq = seq[:-len_ad]
      return seq
    ad = ad[:-1]
  return seq

'''
#NEW FUNCTION
'''
def trim_adapter (seq, ad):
  '''
  example:  adapter ad =                  TGGAATTCTCGGGTGCCAAGGAACTC
            seq =        NTACCGATCTGAGCCATTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACN
            return =     NTACCGATCTGAGCCAT
  '''
  while len(ad) > 6:
    len_ad = len(ad)
    pos = find_str(seq, ad)
    if pos > 0: return seq[:pos]
    else: ad = ad[:-1]
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

def getGenome__ (genome_path, file_ext):
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
  '''
  modified version of original getGenome__()
  this new function can take in either All as the entire genome, or one chromosome at a time.
  this is a measure for reducing memory use
  '''
  genome = dict()
  if chromosomeName == 'All':
    files = [each for each in os.listdir(genome_path) if each.endswith(file_ext)]
  else:
    files = [ chromosomeName + file_ext ]

  for namefile in files :
    file = genome_path + namefile
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
    ## in: ('seq', [freq, nbLoc, ['strd','chr',posChr], ['priSeq',posMirPri,'priFold', 'mkPred','mkStart','mkStop'], ['preSeq',posMirPre,'preFold','mpPred','mpScore'], totalfrq]) 
    fh_out = open (outfile, 'w')

    ## out:
    line = '\t'.join('miRNAseq, frq, nbLoc, strand, chromo, posChr, mkPred, mkStart, mkStop, preSeq, posMirPre, newfbstart, newfbstop, preFold, mpPred, mpScore, totalFrq'.split(', ')) 
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
      data = [miRNAseq, frq, nbLoc, strand, chromo, posChr, mkPred, mkStart, mkStop, preSeq, posMirPre, newfbstart, newfbstop, preFold, mpPred, mpScore, totalFrq] 

      line = '\t'.join([str(d) for d in data])
      print(line, file=fh_out)
    fh_out.close()


def writeTimeLibToFile (timeDict, outfile, appId, paramDict):
  import datetime
  
  totalTimeSec = reduce((lambda x,y : x + y), timeDict.values())
  totalTimeHMS = str(datetime.timedelta(seconds=totalTimeSec))
  totalTimeSec = format(totalTimeSec, '.3f')
  
  fh_out = open (outfile, 'w')
  print(datetime.datetime.now(), file=fh_out)
  print("# Application ID " + appId +"\n", file=fh_out)
  print("Total \t"+totalTimeHMS+ "\t" + totalTimeSec, file=fh_out)

  for lib in timeDict :
    timeLibSec = timeDict[lib]
    timeLibHMS = str(datetime.timedelta(seconds=timeLibSec))
    timeLibSec = format(timeLibSec, '.3f')
    print("Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec, file=fh_out)
    print ("Lib "+lib+"\t"+timeLibHMS+"\t"+timeLibSec) #= stdout
  
  print("\n# SPARK configuration:", file=fh_out)

  for key in paramDict :
    if key.startswith("sc_"):
      print("# " + key + ": " + paramDict[key], file=fh_out)

  print("\n# MirLibSpark configuration:", file=fh_out)

  for key in sorted(paramDict.keys()) :
    if not key.startswith("sc_"):
      print("# " + key + ": " + paramDict[key], file=fh_out)
 
  fh_out.close()


###########################
## WORK IN PROGRESS ==> precursor geno loci summary
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
  print(seqListLine, file=fh_out)
  print(seqListLine, file=fh_out2)

  for k in sorted(dictLibSeqFreq.keys()):
    v = dictLibSeqFreq[k]
    line = k + '\t'
    for i in v: line += str(i) + '\t'
    line = line.rstrip('\t')
    print(line, file=fh_out)

  for k in sorted(dictLibSeqFreq.keys()):
    v = dictLibSeqFreq[k]
    line = k + '\t'
    for i in v:
      if i > 0: i = 1
      line += str(i) + '\t'
    line = line.rstrip('\t')
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

  return sorted(master_distinctPrecursor_infos), sorted(master_predicted_distinctMiRNAs)


def writeTargetsToFile (mirna_and_targets, rep_output, appId):
  ## in:  [miRNAseq, [targets], mirnazipindex]
  ## out: [miRNAseq, t1,t2,t3,t4]
  outfile = rep_output + appId + '_mirna_and_targets.txt'
  fh_out = open (outfile, 'w')
  outfile2 = rep_output + appId + '_mirna_and_topscoredTargets.txt'
  fh_out2 = open (outfile2, 'w')
  topKscored = 5

  for i in mirna_and_targets:
    print(i[0], i[1], str(i[2]).zfill(4), file=fh_out)

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
        targetcollect.append( target + ' ('+ str(score_cur).split('.')[0] +')' )
        seen.append(target)
    #data = [mirnaseq, targetcollect[0], ','.join(targetcollect)]
    data = [mirnaseq, ','.join(targetcollect)]
    line = '\t'.join(data)
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
    print(line, file=fh_out)
  fh_out.close()
  return newDATA 
  
def write_index (data, rep_output, appId):
  #=serial, miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore
  outfile = rep_output + appId + '_precursorindex.txt'
  fh_out = open (outfile, 'w')
  for i in data:
    i[0] = str(i[0]).zfill(4)
    line = '\t'.join( [str(x) for x in i] )
    print(line, file=fh_out)
  fh_out.close()
  __write_html (data, rep_output, appId)

def __write_html (DATA, rep_output, appId):
  #=serial, miRNAseq, strand, chromo, posChr, preSeq, posMirPre, preFold, mkPred, newfbstart, newfbstop, mpPred, mpScore
  infile = rep_output + appId + '_precursorindex.txt'
  outfile = rep_output + appId +'_precursorindex.html'
  fh_out=open(outfile,'w')

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
 
    #path = rep_output + appId + '_' + serial.zfill(4) + '_' + chromo + '_' + poschromo + '.jpg'
    path = appId + '_' + serial.zfill(4) + '_' + chromo + '_' + poschromo + '.jpg'

    l='  <tr>';print(l, file=fh_out)
    l="    <td rowspan=3 style='width: 120px;'><strong>"+ serial + "</strong></td>";print(l, file=fh_out)
    #l="    <td rowspan=3 style='width: 120px;'>"+newid+"</td>";print(l, file=fh_out)
    l="    <td rowspan=3 > <a href="+ path + " target='_blank'><img class='heightSet' src="+ path + " alt='Struture'></a></td>";print(l, file=fh_out)
    l="    <td style='width: 400px;'>"+ mirna + "</td>";print(l, file=fh_out)
    l="    <td> Chr "+ chromo + ":"+ poschromo + " ["+ strand + "] </td>";print(l, file=fh_out)
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
  import string, random
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
  ''' defunct '''
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

def __read_diffguide (infile):
  #= a / b = numerator / denominator
  #= [[numerator, denominator], [numerator, denominator], ...]
  with open (infile, 'r') as fh: diffguide = [x.rstrip('\n').split('->') for x in fh.readlines()]
  needed_infilenames = []
  for i in diffguide[1:]:
    name1 = i[0].split('.')[0]
    name2 = i[1].split('.')[0]
    if name1 not in needed_infilenames: needed_infilenames.append(name1)
    if name2 not in needed_infilenames: needed_infilenames.append(name2)
  return diffguide[1:], needed_infilenames

def __write_diff_output (a, b, rep, appId):
  #rep = '../output/'
  #appId = 'local-1538031347138'
  infile = rep + appId + '_summaryFreq.txt'

  #= a/b: a: numerator, b: denominator
  #a = 'fake_a3'; b = 'fake_a'
  outfile = rep + appId + '_diff_' + a + '_' + b + '.txt'
  with open (infile, 'r') as fh: DATA = [x.rstrip('\n').split('\t') for x in fh.readlines()]

  #= out: Sequence	Iso8S_y2010_2	Iso8S_y2010_1	Fold_change	Z_score	p_value	BH_p_value	Diff_exp
  for i in range( len(DATA[0]) ):
    if DATA[0][i] == a: index_a = i
    if DATA[0][i] == b: index_b = i

  data = [ [x[0], x[index_a], x[index_b]] for x in DATA ]
  title = 'Sequence,' + a + ',' + b + ',Fold_change,Z_score,p_value,BH_p_value,Diff_exp'.split(',')

  import diffAnlysis as dif
  DATA = dif.main (data, title)
  fh_out = open (outfile, 'w')
  for i in DATA: print( '\t'.join([str(x) for x in i]), file=fh_out)
  fh_out.close()


def diff_output (diffguide_file, rep, appId):
  diffguide, _ = __read_diffguide(diffguide_file)
  diff_outs = []
  for i in diffguide:
    a = i[0].split('.')[0]
    b = i[1].split('.')[0]
    __write_diff_output (a, b, rep, appId)
    diff_outs.append( appId + '_diff_' + a + '_' + b )
  return diff_outs


def __dictPathwayDescription (infile):
  #= {ko04978: Mineral absorption	Organismal Systems}, {GO:0006351: transcription, DNA-templated}
  d_pathway_desc = {}
  with open (infile, 'r') as fh: data = [x.rstrip('\n').split('\t') for x in fh.readlines()]
  for i in data: d_pathway_desc[i[0]] = i[1]
  return d_pathway_desc

def __write_namecodefile (folder, ID, diff_outs):
  outfile = folder + '/namecode.txt'
  with open (outfile, 'w') as fh:
    print('background_' + ID + '\tbackground' , file = fh)
    for i in diff_outs:
      base = i.split('_diff_')[1]
      print(i + '\t' + base , file = fh)
  fh.close()

def __create_background (outfile, list_mirna_and_topscoredTargetsKEGGpathway, dict_pathway_description):
  fh_out = open (outfile, 'w')
  for i in list_mirna_and_topscoredTargetsKEGGpathway:
    mirna = i[0]
    pathways = i[1].split(',')
    for p in pathways:
      if p in dict_pathway_description.keys():
        desc = dict_pathway_description[p]
      else: desc = 'to be retrived from KEGG'
      line = '\t'.join( [mirna, p, desc])
      print(line, file=fh_out)
  fh_out.close()

def __create_diff_annotation (rep_output, diff_outs, list_mirna_and_topscoredTargetsKEGGpathway, folder):
  dict_mirna_pathways = {}
  for i in list_mirna_and_topscoredTargetsKEGGpathway:
    mirna = i[0]
    pathways = i[1].split(',')
    dict_mirna_pathways[mirna] = pathways
  
  for infile in diff_outs:
    fh_out = open (folder + '/' + infile, 'w')
    infile = rep_output + infile + '.txt'
    with open (infile, 'r') as fh: data = [x.rstrip('\n').split('\t') for x in fh.readlines()][1:]
    for i in data:
      mirna = i[0]
      UPorDOWN = i[7]
      if UPorDOWN == 'NO': continue
      pathways = dict_mirna_pathways[mirna]
      for p in pathways:
        line = '\t'.join( [mirna, p])
        print(line, file=fh_out)
    fh_out.close()

def __create_inputs_for_enrichment_analysis (diff_outs, pathway_description_file, list_mirna_and_topscoredTargetsKEGGpathway, rep_output, appId):
  '''
  '''
  dict_pathway_description = __dictPathwayDescription (pathway_description_file)

  diffKey = 'UP DOWN'.split()
  ID = appId + '_topscoredTargetsKEGGpathway'
  folder = rep_output + ID
  if not os.path.exists(folder): os.makedirs(folder)
  __write_namecodefile (folder, ID, diff_outs)

  enrich_output = rep_output + ID + '/output_comput_enrich'
  if not os.path.exists(enrich_output): os.makedirs(enrich_output)


  outfile = folder + '/background_' + ID
  __create_background (outfile, list_mirna_and_topscoredTargetsKEGGpathway, dict_pathway_description)
  __create_diff_annotation (rep_output, diff_outs, list_mirna_and_topscoredTargetsKEGGpathway, folder)


def perform_enrichment_analysis (diff_outs, pathway_description_file, list_mirna_and_topscoredTargetsKEGGpathway, rep_output, appId, project_path):
  keyword =  appId + '_topscoredTargetsKEGGpathway'
  __create_inputs_for_enrichment_analysis (diff_outs, pathway_description_file, list_mirna_and_topscoredTargetsKEGGpathway, rep_output, appId)
  #cmd = 'perl compute_enrichment.pl ' + rep_output + keyword + '/namecode.txt ' + rep_output + keyword +'/output_comput_enrich 1'
  cmd = 'perl ' + project_path + '/src/'+ 'compute_enrichment.pl ' + rep_output + keyword + '/namecode.txt ' + rep_output + keyword +'/output_comput_enrich/ 1'
  os.system(cmd)





