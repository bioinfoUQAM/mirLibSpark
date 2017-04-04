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
'''

import os, sys
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
  
  inKvfile = inBasename + '.kv.txt'
  hdfsFile = inBasename + '.hkv.txt'
  
  # Separators
  values_sep = "<>"
  keyval_sep = "::"
  
  # Spark context
  sc = ut.pyspark_configuration("local", "mirLibHadoop", "1g")
  
  # Convert the input file to a Key value file
  ut.convert_seq_freq_file_to_KeyValue(infile, inKvfile, values_sep, keyval_sep)
  
  # Save a local file to HDFS system
  ut.convertTOhadoop(inKvfile, hdfsFile)
  
  # Convert the text file to RDD object
  distFile = sc.textFile(hdfsFile)
  
  # print(distFile.collect())
  
  input_rdd = distFile.flatMap(lambda line: line.split())
  
  # Filtering with DustMasker
  # dmask_obj = mru.prog_dustmasker(values_sep, keyval_sep)
  # dmask_rdd = input_rdd.filter(dmask_obj.dmask_filter_rule)
  
  # Mapping with Bowtie
  bowtie_obj = mru.prog_bowtie(values_sep, keyval_sep, b_index)
  bowtie_rdd = input_rdd.map(bowtie_obj.Bowtie_map_rule)
  
  print bowtie_rdd.collect()
  
