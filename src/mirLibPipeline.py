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

from __future__ import print_function
import sys
import os.path
import time
from os import listdir
#
import utils as ut
import mirLibRules as mru


if __name__ == '__main__' :

  if not len(sys.argv) == 4:
    sys.stderr.write('Three arguments required\nUsage: spark-submit mirLibPipeline.py <path to paramfile> <path to the input directory> <path to the output directory>2>/dev/null\n')
    sys.exit()

  paramfile = sys.argv[1]
  mypath = sys.argv[2]
  rep_output = sys.argv[3]
  
  paramDict = ut.readParam (paramfile)

  # Parameters and cutoffs
  my_sep = paramDict['my_sep']                # Separator
  rep_tmp = paramDict['rep_tmp']              # tmp file folder
  # spark configuration
  appMaster = paramDict['sc_master']             #"local" 
  appName = paramDict['sc_appname']              #"mirLibHadoop"
  mstrMemory = paramDict['sc_mstrmemory']        #"4g"
  execMemory = paramDict['sc_execmemory']        #"4g"
  execNb = paramDict['sc_execnb']                #4
  execCores = paramDict['sc_execcores']          #2
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

  # Spark context
  # (appMaster, appName, masterMemory, execMemory, execNb, execCores):
  sc = ut.pyspark_configuration(appMaster, appName, mstrMemory, execMemory, execNb, execCores)
  sc.addPyFile('utils.py')
  sc.addPyFile('mirLibRules.py')
  
  appId = str(sc.applicationId)
  
  # Object for rule functions

  dmask_obj = mru.prog_dustmasker()
  bowtie_obj = mru.prog_bowtie(b_index)
  prec_obj = mru.extract_precurosrs(genome_path, pri_l_flank, pri_r_flank, pre_flank)
  rnafold_obj = mru.prog_RNAfold()
  mircheck_obj = mru.prog_mirCheck(mcheck_param)
  profile_obj = mru.prog_dominant_profile()

  # fetch library files in mypath
  infiles = [f for f in listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
  
  # time processing of libraries
  timeDict = {}
  
  for infile in infiles :
    print ("--Processing of the library: ", infile)
    
    infile = mypath+infile
    inBasename = os.path.splitext(os.path.basename(infile))[0]
    inKvfile = rep_tmp + inBasename + '.kv.txt'
    hdfsFile = inBasename + '.hkv.txt'

    # Convert the input file to a Key value file
    ut.convert_seq_freq_file_to_KeyValue(infile, inKvfile, my_sep)
  
    # Save a local file to HDFS system
    ut.convertTOhadoop(inKvfile, hdfsFile)
    
    #
    print ("  Start of the processing...", end="\r")
    startLib = time.time()
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

    # Filtering structure with branched loop
    one_loop_rdd = pri_vld_rdd.filter(lambda elem: ut.containsOnlyOneLoop (  elem[1][4][2][ int(elem[1][4][4]) : int(elem[1][4][5])+1 ] ))

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

    #
    endLib = time.time()
    print ("  End of the processing     ", end="\r")
    
    # write results to a file
    outFile = rep_output + inBasename + '_miRNAprediction.txt'
    ut.writeToFile (results, outFile)
    
    timeDict[inBasename] = endLib - startLib
    
  sc.stop() #allow to run multiple SparkContexts
  
  # print executions time  to a file
  outTime = rep_output + appId + '_time.txt'
  ut.writeTimeLibToFile (timeDict, outTime, appId, paramDict)
