'''
program: mirLibPipeline.py
author: M.A.Remita
author: Chao-Jung Wu
date: 2017-03-28
version: 0.00.01

Le programme implemente le pipeline d'analyse des sRAN et prediction des miRNAs. Les principales etapes de predictions sont :
  - 1 Filtrage
  - 2 Alignement 
  - 3 Extraction du precurseur
  - 4 Repliement du precurseur
  - 5 prediction/validation du miRNA
  - 6 validaiton/filtrage avec l'expression
  
La version actuelle accepte un seul argument qui est le fichier contenant les sequences reads a traiter.

personal note - julie's launching command line in her computer:

spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/a_th_3.txt 2>/dev/null

spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/a_th_10.txt 2>/dev/null

$ time spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/110.txt 1>170406_1result_pipeline.txt 2>/dev/null
around 10 mins (6m, 12m, 9.5min)

$ time spark-submit mirLibPipeline.py a_thaliana /home/cloudera/workspace/miRNA_predictor/sudoData/1102.txt 1>/home/cloudera/workspace/miRNA_predictor/logOutput/170407_result_pipeline_1.txt 2>/dev/null
'''

import sys
import os.path
#
import utils as ut
import mirLibRules as mru

if __name__ == '__main__' :

  if not len(sys.argv) == 3:
    sys.stderr.write('Two arguments required\nUsage: spark-submit mirLibPipeline.py <bowtie_index_name> <path to your input> 2>/dev/null\n')
    sys.exit()
  
  b_index = sys.argv[1]
  infile = sys.argv[2]
  
  inBasename = os.path.splitext(os.path.basename(infile))[0]
  
  inKvfile = '../tmp/' + inBasename + '.kv.txt'
  hdfsFile = inBasename + '.hkv.txt'
  
  # Separators
  my_sep = ","

  # Parameters and cutoffs
  #
  genome_path = "../input/ATH/TAIR/Genome/"
  #
  limit_freq = 5 # exclude RNA freq < limit_freq
  limit_len = 14  # exclude RNA length < limit_len
  limit_nbLoc = 5 # exculde nbLoc mapped with bowtie  > limit_nbLoc
  # pri-mirna
  l_flank = 10
  r_flank = 50
  
  # Spark context
  sc = ut.pyspark_configuration("local", "mirLibHadoop", "2g") # yarn-client
  sc.addPyFile('utils.py')
  sc.addPyFile('mirLibRules.py')
  
  # Convert the input file to a Key value file
  ut.convert_seq_freq_file_to_KeyValue(infile, inKvfile, my_sep)
  
  # Save a local file to HDFS system
  ut.convertTOhadoop(inKvfile, hdfsFile)
  
  # Object fo rule functions

  dmask_obj = mru.prog_dustmasker()
  bowtie_obj = mru.prog_bowtie(b_index)
  primir_obj = mru.extract_precurosrs(genome_path, l_flank, r_flank)

  sudo_obj = mru.sudo()
  RNAfold_obj = mru.prog_RNAfold()
  mirCheck_obj = mru.prog_mirCheck()

  
  # Convert the text file to RDD object
  distFile = sc.textFile(hdfsFile)
  input_rdd = distFile.map(lambda line: mru.rearrange_rule(line, my_sep))
  
  # Filtering low frequency
  rm_low_rdd = input_rdd.filter(lambda elem: int(elem[1][1]) > limit_freq)

  # Filtering short length 
  rm_short_rdd = rm_low_rdd.filter(lambda elem: len(str(elem[1][0])) > limit_len)

  # # Filtering with DustMasker
  dmask_rdd = rm_short_rdd.filter(dmask_obj.dmask_filter_rule)

  # Mapping with Bowtie
  bowtie_rdd = rm_low_rdd.map(bowtie_obj.Bowtie_map_rule)

  # Filtering high nbLocations and zero location
  nbLoc_rdd = bowtie_rdd.filter(lambda elem: len(elem[1][2]) > 0 and len(elem[1][2]) < limit_nbLoc)
  
  # Extraction of the pri-miRNA
  # primir_rdd = nbLoc_rdd.map(primir_obj.extract_prec_rule)
  
  #print nbLoc_rdd.collect()

  # sudo_extraction
  sudo_rdd = nbLoc_rdd.map(sudo_obj.sudo_long)
  #print sudo_rdd.collect()

  # RNAfold
  yy_rdd = sudo_rdd.map(RNAfold_obj.RNAfold_map_rule)
  print yy_rdd.collect()

  ## mirCheck
  #zz_rdd = yy_rdd.map(mirCheck_obj.mirCheck_map_rule).filter(mirCheck_obj.mirCheck_filter_rule)
  #print zz_rdd.collect()
