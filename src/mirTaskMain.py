'''
program: mirSparkTask.py
author: Chao-Jung Wu
date: 2017-04-20
version: 0.00.01
'''
import os.path

import utils as ut
import mirLibRules as mru

def SparkTask(paramDict, infile):

  # Parameters and cutoffs
  my_sep = paramDict['my_sep']                # Separator
  rep_tmp = paramDict['rep_tmp']              # tmp file folder
  # spark parameter
  master = paramDict['master']                #"local" 
  appname = paramDict['appname']              #"mirLibHadoop"
  memory = paramDict['memory']                #"2g"
  # genome
  genome_path = paramDict['genome_path']      #"../input/ATH/TAIR/Genome/"
  # cutoffs
  limit_freq = int(paramDict['limit_freq'])   #200      # exclude RNA freq < limit_freq
  limit_len = int(paramDict['limit_len'])     #18       # exclude RNA length < limit_len
  limit_nbLoc = int(paramDict['limit_nbLoc']) #2        # exculde nbLoc mapped with bowtie  > limit_nbLoc
  # bowtie
  b_index = paramDict['b_index']
  # pri-mirna
  pri_l_flank = int(paramDict['pri_l_flank']) #120
  pri_r_flank = int(paramDict['pri_r_flank']) #60
  pre_flank = int(paramDict['pre_flank'])     #30
  # mircheck parameter
  mcheck_param = paramDict['mcheck_param']    #'def'     # def : default parameters / mey : meyers parameters

  inBasename = os.path.splitext(os.path.basename(infile))[0]
  inKvfile = rep_tmp + inBasename + '.kv.txt'
  hdfsFile = inBasename + '.hkv.txt'

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
  profile_obj = mru.prog_dominant_profile()
  
  # Convert the text file to RDD object
  distFile = sc.textFile(hdfsFile)
  input_rdd = distFile.map(lambda line: mru.rearrange_rule(line, my_sep))
  
  # Filtering low frequency
  rm_low_rdd = input_rdd.filter(lambda elem: int(elem[1][1]) > limit_freq)

  # Filtering short length 
  rm_short_rdd = rm_low_rdd.filter(lambda elem: len(str(elem[1][0])) > limit_len)

  # Filtering with DustMasker
  dmask_rdd = rm_short_rdd.filter(dmask_obj.dmask_filter_rule)
  
  # Mapping with Bowtie
  bowtie_rdd = dmask_rdd.map(bowtie_obj.Bowtie_map_rule).persist()

  # Filtering high nbLocations and zero location
  nbLoc_rdd = bowtie_rdd.filter(lambda elem: elem[1][2] > 0 and elem[1][2] < limit_nbLoc)
  
  # Extraction of the pri-miRNA
  primir_rdd = nbLoc_rdd.flatMap(prec_obj.extract_prim_rule)
   
  # pri-miRNA folding
  pri_fold_rdd = primir_rdd.map(lambda elem: rnafold_obj.RNAfold_map_rule(elem, 4))
  
  # Validating pri-mirna with mircheck
  pri_vld_rdd = pri_fold_rdd.map(lambda elem: mircheck_obj.mirCheck_map_rule(elem, 4)).filter(lambda elem: any(elem[1][4]))

  # Filtering structure with branched loop  ####################################################
  one_loop_rdd = pri_vld_rdd.filter(lambda elem: ut.containsOnly1loop (  elem[1][4][2][ int(elem[1][4][4]) : int(elem[1][4][5])+1 ] ))

  # Extraction of the pre-miRNA
  premir_rdd = one_loop_rdd.map(prec_obj.extract_prem_rule)

  # pre-miRNA folding
  pre_fold_rdd = premir_rdd.map(lambda elem: rnafold_obj.RNAfold_map_rule(elem, 5))
  
  # Validating pre-mirna with mircheck
  pre_vld_rdd = pre_fold_rdd.map(lambda elem: mircheck_obj.mirCheck_map_rule(elem, 5)).filter(lambda elem: any(elem[1][5]))

  # you can use chromo_strand as key to search bowtie blocs in the following dict
  dict_bowtie_chromo_strand = ut.return_Bowtie_strandchromo_dict (bowtie_rdd.collect())

  # results of miRNA prediction
  miRNA_rdd = pre_vld_rdd.filter(lambda elem: profile_obj.functionX(elem, dict_bowtie_chromo_strand) )
  results = miRNA_rdd.collect()
  print results
  print len(results)
  sc.stop() #allow to run multiple SparkContexts



