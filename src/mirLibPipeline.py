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
  values_sep = "<>"
  keyval_sep = "::"

  # Cutoffs
  limit_freq = 5 # exclude RNA freq < limit_freq
  limit_len = 14 # exclude RNA length < limit_len
  limit_nbLoc = 5 # exculde nbLoc mapped with bowtie  > limit_nbLoc
  
  # Spark context
  sc = ut.pyspark_configuration("yarn-client", "mirLibHadoop", "2g") # local
  sc.addPyFile('utils.py')
  sc.addPyFile('mirLibRules.py')
  
  # Convert the input file to a Key value file
  ut.convert_seq_freq_file_to_KeyValue(infile, inKvfile, values_sep, keyval_sep)
  
  # Save a local file to HDFS system
  ut.convertTOhadoop(inKvfile, hdfsFile)
  
  # Object fo rule functions
  filter_obj = mru.filter_rules(values_sep, keyval_sep, limit_freq, limit_len, limit_nbLoc)
  dmask_obj = mru.prog_dustmasker(values_sep, keyval_sep)
  bowtie_obj = mru.prog_bowtie(values_sep, keyval_sep, b_index)
  #= insert extraction obj here =#
  sudo_obj = mru.sudo(values_sep, keyval_sep)
  RNAfold_obj = mru.prog_RNAfold(values_sep, keyval_sep)
  mirCheck_obj = mru.prog_mirCheck(values_sep, keyval_sep)
  
  # Convert the text file to RDD object
  distFile = sc.textFile(hdfsFile)
  input_rdd = distFile.flatMap(lambda line: line.split())

  # Filtering low frequency
  rm_low_rdd = input_rdd.filter(filter_obj.lowfreq_filter_rule)

  # Filtering short length
  rm_short_rdd = rm_low_rdd.filter(filter_obj.shortlen_filter_rule)

  # Filtering with DustMasker
  dmask_rdd = rm_short_rdd.filter(dmask_obj.dmask_filter_rule)

  # Mapping with Bowtie
  bowtie_rdd = dmask_rdd.map(bowtie_obj.Bowtie_map_rule)

  # Filtering high nbLocations and zero location
  nbLoc_rdd = bowtie_rdd.filter(filter_obj.nbLocations_filter_rule)

  #print nbLoc_rdd.collect()

  # sudo_extraction
  sudo_rdd = nbLoc_rdd.map(sudo_obj.sudo_long)
  print sudo_rdd.collect()
  # RNAfold
  #yy_rdd = xx_rdd.map(RNAfold_obj.RNAfold_map_rule)

  ## mirCheck
  #zz_rdd = yy_rdd.map(mirCheck_obj.mirCheck_map_rule).filter(mirCheck_obj.mirCheck_filter_rule)



