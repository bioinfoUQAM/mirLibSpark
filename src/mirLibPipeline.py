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

##########################
spark-submit mirLibPipeline.py ../paramfile_julie.txt /home/cloudera/workspace/miRNA_predictor/sudoData/a_th_3.txt 2>/dev/null
'''

import sys
import os.path
#
import utils as ut
import mirLibRules as mru

if __name__ == '__main__' :

  if not len(sys.argv) == 3:
    sys.stderr.write('Two arguments required\nUsage: spark-submit mirLibPipeline.py <path to paramfile> <path to your input> 2>/dev/null\n')
    sys.exit()


  paramfile = sys.argv[1]
  infile = sys.argv[2]
 
  paramDict = ut.readparam (paramfile)
  print paramDict
  # Parameters and cutoffs
  # Separators
  my_sep = paramDict['my_sep']
  # tmp file folder
  rep_tmp = paramDict['rep_tmp']

  # spark parameter
  master = paramDict['master'] #"local" 
  appname = paramDict['appname'] #"mirLibHadoop"
  memory = paramDict['memory'] #"2g"
  # genome
  genome_path = paramDict['genome_path'] #"../input/ATH/TAIR/Genome/"
  # cutoffs
  limit_freq = int(paramDict['limit_freq']) #200            # exclude RNA freq < limit_freq
  limit_len = int(paramDict['limit_len']) #18              # exclude RNA length < limit_len
  limit_nbLoc = int(paramDict['limit_nbLoc']) #2             # exculde nbLoc mapped with bowtie  > limit_nbLoc
  # bowtie
  b_index = paramDict['b_index']
  # pri-mirna
  pri_l_flank = int(paramDict['pri_l_flank']) #120
  pri_r_flank = int(paramDict['pri_r_flank']) #60
  pre_flank = int(paramDict['pre_flank']) #30
  # mircheck parameter
  mcheck_param = paramDict['mcheck_param'] #'def'        # def : default parameters / mey : meyers parameters

  
  inBasename = os.path.splitext(os.path.basename(infile))[0]
  
  inKvfile = rep_tmp + inBasename + '.kv.txt'
  hdfsFile = inBasename + '.hkv.txt'
  print inKvfile
  print hdfsFile

  # Spark context
  sc = ut.pyspark_configuration(master, appname, memory) # yarn-client
  sc.addPyFile('utils.py')
  sc.addPyFile('mirLibRules.py')
  
  # Convert the input file to a Key value file
  ut.convert_seq_freq_file_to_KeyValue(infile, inKvfile, my_sep)
  
  # Save a local file to HDFS system
  ut.convertTOhadoop(inKvfile, hdfsFile)
  
  # Object for rule functions

  dmask_obj = mru.prog_dustmasker()
  bowtie_obj = mru.prog_bowtie(b_index)
  prec_obj = mru.extract_precurosrs(genome_path, pri_l_flank, pri_r_flank, pre_flank)
  rnafold_obj = mru.prog_RNAfold()
  mircheck_obj = mru.prog_mirCheck(mcheck_param)
  
  # Convert the text file to RDD object
  distFile = sc.textFile(hdfsFile)
  input_rdd = distFile.map(lambda line: mru.rearrange_rule(line, my_sep))
  
  # Filtering low frequency
  rm_low_rdd = input_rdd.filter(lambda elem: int(elem[1][1]) > limit_freq)

  # Filtering short length 
  rm_short_rdd = rm_low_rdd.filter(lambda elem: len(str(elem[1][0])) > limit_len)

  # # Filtering with DustMasker
  dmask_rdd = rm_short_rdd.filter(dmask_obj.dmask_filter_rule)
  print dmask_rdd.collect()
  # Mapping with Bowtie
  bowtie_rdd = dmask_rdd.map(bowtie_obj.Bowtie_map_rule)
  

  # Filtering high nbLocations and zero location
  nbLoc_rdd = bowtie_rdd.filter(lambda elem: len(elem[1][2]) > 0 and len(elem[1][2]) < limit_nbLoc)
  
  # Extraction of the pri-miRNA
  primir_rdd = nbLoc_rdd.map(prec_obj.extract_prim_rule)
  
  # pri-miRNA folding
  pri_fold_rdd = primir_rdd.map(lambda elem: rnafold_obj.RNAfold_map_rule(elem, 3))
  
  # Validating pri-mirna with mircheck
  pri_vld_rdd = pri_fold_rdd.map(lambda elem: mircheck_obj.mirCheck_map_rule(elem, 3)).filter(lambda elem: any(elem[1][3]))
  
  # Extraction of the pre-miRNA
  premir_rdd = pri_vld_rdd.map(prec_obj.extract_prem_rule)
  
  # pre-miRNA folding
  pre_fold_rdd = premir_rdd.map(lambda elem: rnafold_obj.RNAfold_map_rule(elem, 4))
  
  # Validating pri-mirna with mircheck
  pre_vld_rdd = pre_fold_rdd.map(lambda elem: mircheck_obj.mirCheck_map_rule(elem, 4)).filter(lambda elem: any(elem[1][4]))
  
  print pre_vld_rdd.collect()

